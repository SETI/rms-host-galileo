##########################################################################################
# host_galileo/oops/ssi.py
##########################################################################################
"""The Galileo SSI host, which makes a Snapshot from each Solid State Imaging image.

The :doc:`user guide </user_guide>` describes the files it reads, how it merges the two
labels of each product, and the errors it raises.
"""

import collections
import re

import julian
import cspyce

from filecache import FCPath
from pdsparser import Pds3Label
from vicar     import VicarError, VicarImage, VicarLabel

import oops
from oops.host._pds3_support import _resolve_pds3_filename
from ._galileo import _Galileo
from .         import _ssi_labels

__all__ = ['GalileoSSI']

# The SSI camera frame is the scan platform frame, whose C kernels are sampled at discrete
# spacecraft clock ticks. SSI images are spaced as closely as 1 unit in the file name,
# which corresponds to 80 clock ticks. Therefore, we use a tolerance of +/-40.
_SPICE_FRAME_NAME = 'GLL_SCAN_PLATFORM'
_SPICE_FRAME_ID = -77001
_TICK_TOLERANCE = 40

# The NAIF ID of the object whose C kernels give the pointing, the scan platform, and of
# the spacecraft clock those C kernels use
_CK_OBJECT_ID = -77001
_CK_SCLK_ID = -77

# The NAIF ID of the SSI instrument, as used in its instrument kernel
_NAIF_ID = -77036

# The PDS3 parser options of this host, which the caller's options update. Every current
# SSI label parses with the strict method, while the fast method misreads the STAR_WINDOW
# sets of 122 labels.
_PDS3_DEFAULTS = {'method': 'strict'}

# The CUT_OUT_WINDOW value of a supplemental index row whose image has no cutout window,
# as 10,791 of the 14,662 current rows give it
_NO_WINDOW = [-1, -1, -1, -1]

# The dimensions of an SSI image in pixels, as the instrument kernel must give them
_SSI_PIXELS = 800

# Seconds from the start of the frame named by SPACECRAFT_CLOCK_START_COUNT to the opening
# of the shutter, keyed by FRAME_DURATION in seconds. These are the median offsets between
# the start of the exposure (IMAGE_TIME, the middle of the exposure, minus half of
# EXPOSURE_DURATION) and the clock count across the 14,585 images of GO_0002-GO_0023 that
# give both; about 90% of images agree to within 1/3 s, the resolution of a count given to
# the mod-91 field.
_SHUTTER_DELAY = {
     2.333: 0.225,
     8.667: 1.149,
    15.167: 1.150,
    30.333: 2.817,
    60.667: 6.483,
}
_FRAME_DURATION_TOLERANCE = 0.05

# The RIM count of the first frame in the archive taken under the phase-2 flight software
# of the Jupiter tour, GO_0017 J0/OPNAV/C0346405900R. Every earlier frame carries the
# phase-1 FIBE flags in its VICAR header, the last of them at RIM 2492218, and every
# later frame the phase-2 MOFIBE flags. Under phase 1, an EXTENDED exposure opens its
# shutter one frame earlier, in the frame before the one the count names: across the
# 14,672 frames of GO_0002-GO_0023, the 160 phase-1 EXTENDED frames of 30.333 s that give
# IMAGE_TIME open at a median of -27.184 s, the normal 2.817 s less one frame. A phase-2
# EXTENDED exposure opens at about the normal delay: a median of 2.816 s for the 81 such
# frames of 30.333 s. All 77 frames with no IMAGE_TIME are phase-1 frames, so the rule
# changes no current product, but it keeps tstart_from_sclk() right for every frame.
_PHASE2_FIRST_RIM = 3464059

# The separators between the fields of a clock count, in the label form "RRRRRRRR.MM" or
# the SPICE form "RRRRRRRR:MM:S:T", as _Galileo.tdb_from_sclk() splits it
_SCLK_SEPARATORS = re.compile(r'[^0-9A-Za-z]+')

# The subfields that give a metadata item unchanged, and the items' names. Each subfield
# is inserted only if the item has a value.
_ITEM_SUBFIELDS = [
    ('filter'                      , 'FILTER_NAME'                 ),
    ('label_target'                , 'TARGET_NAME'                 ),
    ('observation_id'              , 'OBSERVATION_ID'              ),
    ('spacecraft_clock_start_count', 'SPACECRAFT_CLOCK_START_COUNT'),
]

# Target name as found in PDS label -> correct name. 150 labels spell BLACK_SKY and
# NON_SCIENCE with a space.
_TARGET_NAME_REPAIRS = {
    'BLACK_SKY'  : 'NONE',
    'BLACK SKY'  : 'NONE',
    'NON_SCIENCE': 'NONE',
    'NON SCIENCE': 'NONE',
    'STAR'       : 'NONE',
    'HD151288'   : 'HD 151288',
    'J RINGS'    : 'JUPITER_RING_PLANE',
}

# These are star targets and need to be defined as LightSources
# "BET_ARI", "DEL_PSC", "GAM_ORI", "HD 151288", "PLEIADES", "SIRIUS", "TAU_CET"
_TARGET_STARS = re.compile(r'(HD 151288|PLEIADES|SIRIUS|[A-Z]{3}_[A-Z]{3})$')

# What the files of one product give: the merged metadata items, the Snapshot inputs
# taken from them, the two parsed labels as dictionaries (each None if not read), the
# image array (None if not read), and the paths of the image file and of the label (None
# if not read)
_Product = collections.namedtuple('_Product', ['items', 'inputs', 'label_dict',
                                               'header_dict', 'data', 'image_path',
                                               'label_path'])

# The Snapshot inputs taken from the metadata: the exposure in seconds; the midtime in
# seconds TDB, or None if it comes from the clock; the clock count and the seconds from
# the start of its frame to the opening of the shutter, or None if the midtime is known;
# the FOV of the telemetry format; the cutout window as zero-based [l0, s0, l1, s1], or
# None; the target; and the light source
_Inputs = collections.namedtuple('_Inputs', ['texp', 'tmid', 'sclk', 'fov', 'window',
                                             'target', 'lightsource'])


class GalileoSSI(oops.Host):
    """The Galileo SSI host: Snapshot constructors for Solid State Imaging images."""

    NAME = 'Galileo SSI'

    # A survey of all 23 SSI metadata indexes (GO_0002-GO_0023 and the cumulative
    # GO_0999, 29,324 rows), reading VOLUME_ID and FILE_SPECIFICATION_NAME from both
    # tables, found the same files in the same order in every index and its supplement,
    # none listed twice.

    #: The rows of an SSI index and of its supplement are paired by position, because
    #: every SSI metadata index lists its files in the same order as its supplement.
    MERGE_SUPPLEMENT_BY_FILE = False

    _INSTRUMENT_KERNEL = None
    _FOVS = {}
    _initialized = False

    @staticmethod
    def from_file(fileinfo, *, full_fov=False, return_all_planets=False,
                  pds3_options=None, astrometry=False, select=None, target=None,
                  lightsource=None, path=None, frame=None, cadence=None, fov=None,
                  calibrations=None, timeshift=None, navigation=None, tracker=None,
                  parallel=None):
        """A Snapshot object based on a given Galileo SSI image file.

        Given the VICAR image file, only that file is read. Given its PDS3 label, the
        image file that the label's IMAGE pointer names, in the label's directory and in
        any case, is read as well, unless `astrometry` is True. The array always comes
        from the VICAR file's own header, never from the record number in the label's
        IMAGE pointer, and one read of the image file gives both its header and its
        array.

        The metadata come from every label that was read, the PDS3 label and the VICAR
        header. A header value that the host converts into a PDS3 item, such as EXP into
        EXPOSURE_DURATION, is given only under its PDS3 name, with the type the PDS3
        label gives it. An item that no label read gives is left out, and so is one whose
        value is a PDS3 null text (UNK, N/A, NULL, UNKNOWN, NOT APPLICABLE, "---", or an
        empty string, in any case), a header value of -999, or a label ORBIT_NUMBER
        outside 0-35, the orbits of the Jupiter tour. Every other value is given as the
        label writes it. Where both labels give an item and the values differ, each item
        takes its value from one fixed label for every image: the VICAR header for
        EXPOSURE_DURATION, for the pointing and viewing geometry that the PDS3 label
        rounds or miscomputes (such as DECLINATION, TARGET_CENTER_DISTANCE, and
        SUB_SOLAR_LONGITUDE), and for the compression ratios; the PDS3 label for every
        other item, including CENTRAL_BODY_DISTANCE, COMPRESSION_TYPE, and TARGET_NAME.
        If the winning label lacks a value, the other label's value is used.

        The exposure is centered on IMAGE_TIME. If no label read gives it, the exposure
        starts at the time :meth:`tstart_from_sclk` gives.

        The C kernels are not checked: if they lack the pointing of an image, the error
        of SPICE is raised when the frame is evaluated.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | Pds3Label | VicarLabel | VicarImage):
                The path to the PDS3 label or the VICAR image file, or its parsed content
                as a Pds3Label, VicarLabel, or VicarImage. A path whose name ends in
                ".lbl", in any case, is a PDS3 label; any other path is a VICAR image
                file.
            full_fov (bool, optional): If True, the full 800x800 image and its FOV are
                returned, and the `window` subfield gives the bounds of the cutout
                window; otherwise, the image is trimmed to the cutout window and the FOV
                is a SliceFOV of that window.
            return_all_planets (bool, optional): Include kernels for all planets, not just
                the target of the mission phase.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser. The
                method is "strict" unless they give another.
            astrometry (bool, optional):
                True to specify that the returned Observations will only be used for
                timing and astrometry. In this case, no data arrays are read and returned
                in the Observation.
            select (int | slice | str | tuple[int | slice | str, ...], optional): A
                choice from the image's one Snapshot, made by ``Host._apply_select`` of
                oops as from a tuple holding that Snapshot: 0 and -1 give the Snapshot,
                and a slice gives a tuple. None, the default, gives the Snapshot.
            target (str | oops.Body, optional): Override for the default target Body of
                the observation. Use "NONE" for inertial pointing, indicating that no
                Solar System body was tracked.
            lightsource (oops.LightSource, optional): Override for the default LightSource
                of the observation.
            path (str | oops.Path, optional): Override for the Path of the observer.
            frame (str | oops.Frame, optional): Override for the Frame of the observing
                instrument, used as given.
            cadence (oops.Cadence, optional): Override for the Cadence of the observation.
            fov (oops.FOV, optional): Override for the default FOV of the observing
                instrument.
            calibrations (oops.Calibration | list[oops.Calibration], optional): Override
                for the calibration or list of calibrations.
            timeshift (tuple[float, str], optional): Assign a Fittable time shift to one
                or more attributes of the Observation. The first item is the initial time
                shift in seconds. The second item is a string indicating which attributes
                are to be time-shifted: "T" for the start/stop time, "P" for the path,
                and/or "F" for the frame. If multiple letters are given, the time shifts
                are coupled to whichever attribute is listed first.
            navigation (tuple[float, ...], optional): Two or three rotation angles in
                radians to apply to the frame in order to align the calculated geometry
                with the observed geometry. In this case, a fittable Navigation frame
                wraps the default frame. Specify three angles to include a rotation about
                the optic axis or two for a pointing offset without rotation. Use (0,0) or
                (0,0,0) as the input if you have no starting guess.
            tracker (str, optional): Assign a TrackerFrame to the Observation, to ensure
                that the target body remains at a fixed position within the FOV. Use one
                of "start", "midtime", and "end", indicating the time within the
                Observation at which the target body's calculated position within the FOV
                is accurate.
            parallel (Observation, optional): A parallel Observation (same origin and
                time, different frame and FOV) relative to which this Observation has a
                fixed offset. Note that no more than one of `parallel` and `tracker` can
                be specified.

        Returns:
            Snapshot | tuple[Snapshot, ...]: The observation, or the choice that `select`
            makes, with `data` inserted unless `astrometry` is True, and these subfields:

            - `dict`: the merged metadata described above. A VICAR keyword that the host
              converts appears only under its PDS3 name: RIM and MOD91 as
              SPACECRAFT_CLOCK_START_COUNT ("RRRRRRRR.MM"), the SCET and ERT fields as
              IMAGE_TIME and EARTH_RECEIVED_TIME (each a datetime), RATE as
              FRAME_DURATION, and GAIN as GAIN_MODE_ID; FILTER_NAME comes from FILTER,
              and EXPOSURE_TYPE from the MOFIBE or FIBE flags. MOD10 and MOD8 keep their
              own names, so an astrometry read of the PDS3 label leaves them out.
              CUT_OUT_WINDOW and TRUTH_WINDOW are lists in the label's order.
            - `label_dict`: the PDS3 label as ``Pds3Label.as_dict()`` gives it, if the
              PDS3 label was read.
            - `header_dict`: the VICAR header as ``VicarLabel.as_dict()`` gives it, if the
              image file was read. Its HOST and BHOST are those of the machine that reads
              it, because the VICAR reader replaces the file's values.
            - `host` ("Galileo"), `instrument` ("SSI"), `time_from_sclk` (True if no
              IMAGE_TIME was read, so that the time comes from the spacecraft clock
              count), `window` (the cutout window as zero-based line and sample limits
              ``[l0, s0, l1, s1]``, or the full 800x800 frame if there is none), and
              `full_fov`.
            - `filter` (FILTER_NAME), `label_target` (TARGET_NAME), `observation_id`
              (OBSERVATION_ID), and `spacecraft_clock_start_count`
              (SPACECRAFT_CLOCK_START_COUNT), each only if `dict` gives the item.
            - `spice_kernels`; `spice_to_frame` (the identity), `spice_frame_name`
              ("GLL_SCAN_PLATFORM"), and `spice_frame_id` (-77001); `ck_object_id`
              (-77001) and `ck_sclk_id` (-77), the NAIF IDs of the C kernel object and
              of its clock; and `spice_frame_is_frozen` (False).
            - `filepath` and `basename`, which name the image file also when `fileinfo`
              is its PDS3 label, and `image_url`, its absolute path or URL;
              `image_abspath`, its local path, if the image file was read; and
              `label_url` and `label_abspath` if the PDS3 label was read.

        Raises:
            TypeError: If a keyword argument is not one of those listed.
            oops.host.PdsHostError: If the PDS3 label cannot be parsed with
                `pds3_options`, including a value that the parser cannot convert, or if
                it has no IMAGE pointer. The message names the label.
            oops.host.VicarHostError: If the VICAR image file cannot be read, including
                one shorter than its header says, or if its array has more than one band.
                The message names the file.
            oops.host.HostFormatError: If the array does not have the dimensions that the
                PDS3 label gives it, or if the metadata cannot be interpreted: for
                example, if no label read gives the exposure duration, the telemetry
                format, or a time (IMAGE_TIME, or the spacecraft clock count and frame
                duration), if the telemetry format is not one of SSI's, or if the time
                must come from the spacecraft clock count and FRAME_DURATION is not one
                of the SSI frame durations. The message names the label, or the image
                file if no label was read.
            oops.host.WrongHostError: If the PDS3 label or VICAR image file is not that of
                a Galileo SSI image.
            FileNotFoundError: If a file does not exist, including the image file that
                the PDS3 label names, unless `astrometry` is True.
            OSError: If a file cannot be read.
        """

        GalileoSSI.initialize()  # define everything the first time through; use defaults
                                 # unless initialize() was called explicitly.

        product = GalileoSSI._read_product(fileinfo, astrometry=astrometry,
                                           pds3_options=pds3_options)

        overrides = {
            'target'      : target,
            'lightsource' : lightsource,
            'path'        : path,
            'frame'       : frame,
            'cadence'     : cadence,
            'fov'         : fov,
            'calibrations': calibrations,
            'timeshift'   : timeshift,
            'navigation'  : navigation,
            'tracker'     : tracker,
            'parallel'    : parallel,
        }

        obs = GalileoSSI._make_snapshot(product.items, product.inputs, data=product.data,
                                        full_fov=full_fov,
                                        return_all_planets=return_all_planets,
                                        overrides=overrides)
        oops.Host._set_path_attributes(obs, product.image_path, product.label_path)
        if product.label_dict is not None:
            obs.insert_subfield('label_dict', product.label_dict)
        if product.header_dict is not None:
            obs.insert_subfield('header_dict', product.header_dict)
        return GalileoSSI._select(obs, select)

    @staticmethod
    def _read_product(fileinfo, *, astrometry, pds3_options):
        """Read the files of one SSI product and interpret their metadata.

        Every read of a label or an image file by :meth:`from_file` goes through this
        method, which turns a failure to read or interpret either file into a host error
        that names the file.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | Pds3Label | VicarLabel | VicarImage):
                The argument of :meth:`from_file`.
            astrometry (bool): True to leave the image file of a PDS3 label unread and
                the array of any image unread.
            pds3_options (dict[str, Any] | None): The caller's PDS3 parser options.

        Returns:
            _Product: What the files give.

        Raises:
            oops.host.PdsHostError: If the PDS3 label cannot be parsed or has no IMAGE
                pointer.
            oops.host.VicarHostError: If the VICAR image file cannot be read or its array
                has more than one band.
            oops.host.HostFormatError: If the array does not have the dimensions that the
                PDS3 label gives it, or if the metadata cannot be interpreted.
            oops.host.WrongHostError: If a file is not that of a Galileo SSI image.
            OSError: If a file does not exist or cannot be read.
        """

        label = None
        header = None
        image = None
        if isinstance(fileinfo, Pds3Label):
            label = fileinfo
        elif isinstance(fileinfo, VicarLabel):
            header = fileinfo
        elif isinstance(fileinfo, VicarImage):
            image = fileinfo
        elif FCPath(fileinfo).suffix.lower() == '.lbl':
            label = GalileoSSI._read_label(FCPath(fileinfo), pds3_options)
        else:
            parsed = GalileoSSI._read_vicar(FCPath(fileinfo), astrometry=astrometry)
            if not GalileoSSI._detect_in_vicar(parsed):
                raise oops.Host._no_host_error([GalileoSSI], 'VICAR file', fileinfo)
            if astrometry:
                header = parsed
            else:
                image = parsed

        label_dict = None
        label_path = None
        if label is not None:
            label_dict = label.as_dict()
            label_path = FCPath(label.filepath)
            image_path = GalileoSSI._image_path(label_dict, label_path,
                                                find=not astrometry)
            if not astrometry:
                image = GalileoSSI._read_vicar(image_path, astrometry=False)
        else:
            image_path = FCPath((image if header is None else header).filepath)
            if header is not None and not astrometry:
                image = GalileoSSI._read_vicar(image_path, astrometry=False)

        data = None
        if image is not None:
            header = image.label
            if not astrometry:
                data = GalileoSSI._array_2d(image, image_path)

        header_dict = None if header is None else header.as_dict()

        # Interpreting the metadata is plain Python with no SPICE call, so a KeyError,
        # TypeError, or ValueError here comes from the content of the files
        named = image_path if label_path is None else label_path
        try:
            items = _ssi_labels.merge_items(
                None if label_dict is None else _ssi_labels.label_items(label_dict),
                None if header is None else _ssi_labels.header_items(header))
            inputs = GalileoSSI._interpret(items)
            shape = None
            if data is not None and label_dict is not None:
                image_object = label_dict['IMAGE']
                shape = (image_object['LINES'], image_object['LINE_SAMPLES'])
        except (KeyError, TypeError, ValueError) as err:
            raise oops.host.HostFormatError(f'cannot interpret the Galileo SSI product '
                                            f'{named}: {err}') from err

        # An array of other dimensions than its label gives would silently give the
        # wrong cutout window
        if shape is not None and data.shape != shape:
            raise oops.host.HostFormatError(f'Galileo SSI image file {image_path} holds '
                                            f'a {data.shape[0]}x{data.shape[1]} array, '
                                            f'where its label {label_path} gives '
                                            f'{shape[0]}x{shape[1]}')

        return _Product(items, inputs, label_dict, header_dict, data, image_path,
                        label_path)

    @staticmethod
    def _read_label(filepath, pds3_options):
        """The parsed PDS3 label of an SSI image.

        Parameters:
            filepath (FCPath): The path to the label.
            pds3_options (dict[str, Any] | None): The caller's PDS3 parser options.

        Returns:
            Pds3Label: The label.

        Raises:
            oops.host.PdsHostError: If the label cannot be parsed, including a value that
                the parser cannot convert.
            oops.host.WrongHostError: If the label is not that of a Galileo SSI image.
            OSError: If the label does not exist or cannot be read.
        """

        options = GalileoSSI._pds3_options(pds3_options)
        try:
            return GalileoSSI._read_fileinfo(filepath, formats='P', pds3_options=options)
        except oops.OopsValueError as err:
            raise oops.host.PdsHostError(f'cannot parse the PDS3 label {filepath}: '
                                         f'{err}') from err

    @staticmethod
    def _image_path(label_dict, label_path, *, find):
        """The image file that a PDS3 label describes.

        Parameters:
            label_dict (dict): The parsed label.
            label_path (FCPath): The path to the label.
            find (bool): True to look for the file, which then need not have the same
                case as in the label; False to name it without touching the file system.

        Returns:
            FCPath: The file named by the label's IMAGE pointer, in the directory of the
            label.

        Raises:
            oops.host.PdsHostError: If the label has no IMAGE pointer.
            FileNotFoundError: If `find` is True and the directory holds no file of that
                name in any case.
        """

        if '^IMAGE' not in label_dict:
            raise oops.host.PdsHostError(f'no ^IMAGE pointer in Galileo SSI PDS3 label '
                                         f'{label_path}')
        name = label_dict['^IMAGE']
        if find:
            return _resolve_pds3_filename(label_path.parent, name)
        return label_path.parent / name

    @staticmethod
    def _read_vicar(filepath, *, astrometry):
        """The VICAR header, or the whole VICAR image, of an SSI image file.

        The header is read with the loose VICAR reader: 23 of the 719 GO_0002 RAW_CAL
        headers, such as that of C0003061900R, carry a byte in BARC that is not ASCII,
        which the strict reader refuses.

        Parameters:
            filepath (FCPath): The image file.
            astrometry (bool): True to read only the header; False to read the whole
                image.

        Returns:
            VicarLabel | VicarImage: The header if `astrometry` is True; otherwise, the
            image.

        Raises:
            oops.host.VicarHostError: If the file is not a valid VICAR file, including an
                image shorter than its header says. The message names the file.
            OSError: If the file does not exist or cannot be read.
        """

        # A truncated array fails to take the header's shape with a bare ValueError
        try:
            if astrometry:
                return VicarLabel(filepath, strict=False)
            return VicarImage(filepath, strict=False)
        except (VicarError, ValueError) as err:
            raise oops.host.VicarHostError(f'cannot read the Galileo SSI image file '
                                           f'{filepath}: {err}') from err

    @staticmethod
    def _array_2d(image, filepath):
        """The array of a VICAR image, in two dimensions.

        Parameters:
            image (VicarImage): The image.
            filepath (FCPath): The image file, for the message.

        Returns:
            numpy.ndarray: The array.

        Raises:
            oops.host.VicarHostError: If the array has more than one band.
        """

        try:
            return image.data_2d
        except VicarError as err:
            raise oops.host.VicarHostError(f'cannot read the Galileo SSI image file '
                                           f'{filepath}: {err}') from err

    @staticmethod
    def from_index(filepath, *, supplement=oops.Host.SUPPLEMENT, label=None,
                   row_dicts=None, pds3_options=None, filter=None, full_fov=False,
                   return_all_planets=False, select=None, target=None, lightsource=None,
                   path=None, frame=None, cadence=None, fov=None, calibrations=None,
                   timeshift=None, navigation=None, tracker=None, parallel=None):
        """An iterator over Snapshots, one for each row of a Galileo SSI index file.

        It reads a Ring-Moon Systems Node metadata index together with its supplement,
        which gives TELEMETRY_FORMAT_ID and CUT_OUT_WINDOW. The call itself reads the
        labels of the two tables, so that a wrong call fails at once; the rows are then
        read as the iteration proceeds, 100 at a time, and each row's Snapshot is built
        when the iteration reaches that row.

        Every SSI index lists its files in the same order as its supplement, so each
        index row is paired with the supplement row in the same position, and the two
        rows must give the same file. No other column is compared, because the two
        tables write SPACECRAFT_CLOCK_START_COUNT in different forms. Where both rows
        give a column, the index row's value is kept.

        The C kernels are not checked: if they lack the pointing of an image, the error
        of SPICE is raised when its frame is evaluated.

        Parameters:
            filepath (str | pathlib.Path | FCPath): The full path to the PDS label of a
                Galileo SSI index, such as "GO_0017_index.lbl".
            supplement (str | pathlib.Path | FCPath | object, optional): The path to the
                PDS label of the supplemental index, in any directory. The default,
                `oops.Host.SUPPLEMENT`, reads the file in the directory of `filepath`
                whose name matches except that "_index" is replaced by
                "_supplemental_index" in upper or lower case. The supplement must exist.
            label (Pds3Label, optional): The PDS label of `filepath` if it has already
                been read.
            row_dicts (collections.abc.Iterable[dict[str, Any]], optional): The rows of
                the index, each already merged with its supplement row. If they are given,
                the tables are not read and `supplement` is not used.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser for the
                label of the index; the label of the supplement is parsed with their
                method alone. The method is "strict" unless they give another.
            filter (dict[str, Any], optional): A dictionary of index column names and
                their values or a set or tuple of values. If provided, only rows of the
                index in which each named column equals the value, falls in the set of
                values, or falls in a range defined by a tuple of two values, bounds
                included, are yielded.
            full_fov (bool, optional): If True, each field of view covers the full image
                rather than the cutout window.
            return_all_planets (bool, optional): Include kernels for all planets, not just
                the target of the mission phase.
            select (int | slice | str | tuple[int | slice | str, ...], optional): A
                choice from each row's one Snapshot, as described for :meth:`from_file`.
            target (str | oops.Body, optional): Override for the default target Body of
                every observation. Use "NONE" for inertial pointing, indicating that no
                Solar System body was tracked.
            lightsource (oops.LightSource, optional): Override for the default LightSource
                of every observation.
            path (str | oops.Path, optional): Override for the Path of the observer.
            frame (str | oops.Frame, optional): Override for the Frame of the observing
                instrument, used as given.
            cadence (oops.Cadence, optional): Override for the Cadence of the observation.
            fov (oops.FOV, optional): Override for the FOV of every observation,
                replacing the FOV that each row's telemetry format and cutout window would
                otherwise define.
            calibrations (oops.Calibration | list[oops.Calibration], optional): Override
                for the calibration or list of calibrations.
            timeshift (tuple[float, str], optional): Assign a Fittable time shift to one
                or more attributes of the Observation. The first item is the initial time
                shift in seconds. The second item is a string indicating which attributes
                are to be time-shifted: "T" for the start/stop time, "P" for the path,
                and/or "F" for the frame. If multiple letters are given, the time shifts
                are coupled to whichever attribute is listed first.
            navigation (tuple[float, ...], optional): Two or three rotation angles in
                radians to apply to the frame in order to align the calculated geometry
                with the observed geometry. In this case, a fittable Navigation frame
                wraps the default frame. Specify three angles to include a rotation about
                the optic axis or two for a pointing offset without rotation. Use (0,0) or
                (0,0,0) as the input if you have no starting guess.
            tracker (str, optional): Assign a TrackerFrame to the Observation, to ensure
                that the target body remains at a fixed position within the FOV. Use one
                of "start", "midtime", and "end", indicating the time within the
                Observation at which the target body's calculated position within the FOV
                is accurate.
            parallel (Observation, optional): A parallel Observation (same origin and
                time, different frame and FOV) relative to which this Observation has a
                fixed offset. Note that no more than one of `parallel` and `tracker` can
                be specified.

        Returns:
            collections.abc.Iterator[Snapshot | tuple[Snapshot, ...]]: One observation
            for every row of the index that passes `filter`, in index order, including
            any row with a zero exposure, which has a zero `texp`. A caller that needs
            ``len()`` or indexing wraps the call in ``list()``. Each observation has the
            subfields that :meth:`from_file` describes, except `data`, `label_dict`,
            `header_dict`, and the subfields that name files other than `filepath` and
            `basename`. `dict` is the row, merged with its supplement row, with every
            value as the table reader gives it, including its "_mask" entries.
            `filepath` is the row's VOLUME_ID and FILE_SPECIFICATION_NAME joined by a
            slash, as in "GO_0017/G1/GANYMEDE/C0349632100R.LBL", and `basename` is its
            last component. The exposure is centered on the row's IMAGE_TIME, or starts
            at the time :meth:`tstart_from_sclk` gives if IMAGE_TIME is a PDS3 null text
            such as UNK. A CUT_OUT_WINDOW of -1, -1, -1, -1 means that the image has no
            cutout window.

        Raises:
            TypeError: At call time, if a keyword argument is not one of those listed.
            oops.host.HostIndexError: At call time, if `supplement` is None, if
                `supplement` is `oops.Host.SUPPLEMENT` and the name of `filepath`
                contains no "_index", or if either label does not describe a table that
                can be read. When the iteration reaches a row, if the index row and the
                supplement row in the same position give different files, or if the
                supplement ends before the index; after the last index row, if the
                supplement has rows left over. The table reader raises it for a table
                file that does not hold the rows its label describes.
            oops.host.HostFormatError: When the iteration reaches a row whose time must
                come from the spacecraft clock count while its FRAME_DURATION is not one
                of the SSI frame durations.
            oops.host.PdsHostError: At call time, if the label of the index or of the
                supplement cannot be parsed.
            oops.host.WrongHostError: At call time, if `filepath` is not a Galileo SSI
                index.
            FileNotFoundError: At call time, if the label of the index or of the
                supplement does not exist.
            OSError: At call time, if a label cannot be read; when the iteration starts,
                if a table file does not exist or cannot be read.
        """

        if supplement is None:
            raise oops.host.HostIndexError(f'Galileo SSI index needs its supplemental '
                                           f'index: {filepath}')

        GalileoSSI.initialize()

        # The two tables also share SPACECRAFT_CLOCK_START_COUNT, but the supplement
        # writes it in the SPICE form "RRRRRRRR:MM:S:T" and the index in the label form
        # "RRRRRRRR.MM", so a comparison of every shared column would reject every real
        # pair. No column is compared; the pairing still compares the files of the two
        # rows.
        rows = GalileoSSI._index_rows(filepath, supplement=supplement, label=label,
                                      row_dicts=row_dicts,
                                      pds3_options=GalileoSSI._pds3_options(pds3_options),
                                      filter=filter, validate=())

        overrides = {
            'target'      : target,
            'lightsource' : lightsource,
            'path'        : path,
            'frame'       : frame,
            'cadence'     : cadence,
            'fov'         : fov,
            'calibrations': calibrations,
            'timeshift'   : timeshift,
            'navigation'  : navigation,
            'tracker'     : tracker,
            'parallel'    : parallel,
        }

        return (GalileoSSI._row_snapshot(row, full_fov=full_fov,
                                         return_all_planets=return_all_planets,
                                         select=select, overrides=overrides)
                for row in rows)

    @staticmethod
    def _row_snapshot(row, *, full_fov, return_all_planets, select, overrides):
        """The Snapshot of one row of an index, or the choice that `select` makes.

        The row is used as the table gives it, with no check.

        Parameters:
            row (dict): The row, merged with its supplement row.
            full_fov (bool): True to keep the full FOV; False to trim it to the cutout
                window.
            return_all_planets (bool): Include kernels for all planets, not just the
                target of the mission phase.
            select (int | slice | str | tuple[int | slice | str, ...] | None): The choice
                described for :meth:`from_file`.
            overrides (dict[str, Any]): The values of the standard overrides.

        Returns:
            Snapshot | tuple[Snapshot, ...]: The observation, or the choice from it.
        """

        obs = GalileoSSI._make_snapshot(row, GalileoSSI._interpret(row), data=None,
                                        full_fov=full_fov,
                                        return_all_planets=return_all_planets,
                                        overrides=overrides)
        oops.Host._set_path_attributes(obs, f'{row["VOLUME_ID"]}/'
                                            f'{row["FILE_SPECIFICATION_NAME"]}', None)
        return GalileoSSI._select(obs, select)

    @staticmethod
    def _interpret(dict_):
        """The Snapshot inputs that the metadata of an image or an index row give.

        Nothing here calls SPICE, so the clock count is converted later.

        Parameters:
            dict_ (dict): The merged metadata of an image, or an index row merged with
                its supplement row.

        Returns:
            _Inputs: The inputs.

        Raises:
            oops.host.HostFormatError: If the time must come from the spacecraft clock
                count and FRAME_DURATION is not one of the SSI frame durations.
            KeyError: If an item that is needed is missing, or the telemetry format is
                not one of SSI's.
            TypeError: If a value has the wrong type.
            ValueError: If a value cannot be interpreted.
        """

        # PDS3 labels and index rows use the same names
        texp = dict_['EXPOSURE_DURATION'] / 1000.

        # IMAGE_TIME is the middle of the exposure. Some RAW_CAL frames in GO_0002 and
        # GO_0003 have no IMAGE_TIME, but every label carries the spacecraft clock count.
        # The merged metadata of an image leave a null IMAGE_TIME out, while an index row
        # gives it as the table writes it.
        image_time = dict_.get('IMAGE_TIME')
        if image_time is None or oops.host.is_pds3_null(image_time):
            tmid = None
            sclk = GalileoSSI._sclk_start(dict_)
        else:
            tmid = julian.tdb_from_tai(julian.tai_from_iso(image_time))
            sclk = None

        fov = GalileoSSI._FOVS[dict_['TELEMETRY_FORMAT_ID']]

        # The cutout window is (line, sample, lines, samples), with a one-based origin
        window = dict_.get('CUT_OUT_WINDOW')
        if window is not None:
            window = [int(value) for value in window]
            if window == _NO_WINDOW:
                window = None
            else:
                (l1, s1, nl, ns) = window
                window = [l1 - 1, s1 - 1, l1 + nl - 1, s1 + ns - 1]

        label_target = dict_.get('TARGET_NAME')
        target = _TARGET_NAME_REPAIRS.get(label_target, label_target)
        lightsource = 'SUN'
        if target is not None and _TARGET_STARS.match(target):
            lightsource = oops.lightsource.star_lookup(target)
            target = 'NONE'

        return _Inputs(texp, tmid, sclk, fov, window, target, lightsource)

    @staticmethod
    def _make_snapshot(dict_, inputs, *, data, full_fov, return_all_planets, overrides):
        """A Snapshot built from the metadata of an image or a row of an index.

        The subfields that name files are left to the caller.

        Parameters:
            dict_ (dict): The merged metadata of an image, or an index row merged with
                its supplement row.
            inputs (_Inputs): The inputs that :meth:`_interpret` takes from `dict_`.
            data (numpy.ndarray | None): The full image array; None to omit the `data`
                subfield.
            full_fov (bool): True to keep the full FOV and image; False to trim both to
                the cutout window.
            return_all_planets (bool): Include kernels for all planets, not just the
                target of the mission phase.
            overrides (dict[str, Any]): The values of the standard overrides, `target`,
                `lightsource`, `path`, `frame`, `cadence`, `fov`, `calibrations`,
                `timeshift`, `navigation`, `tracker` and `parallel`, as described for
                :meth:`from_file`.

        Returns:
            Snapshot: The observation.
        """

        if inputs.tmid is None:
            (count, delay) = inputs.sclk
            tstart = _Galileo.tdb_from_sclk(count) + delay
        else:
            tstart = inputs.tmid - inputs.texp / 2.

        # Make sure the SPICE frame is defined unless a custom frame was given
        using_cks = overrides['frame'] is None
        if using_cks:
            GalileoSSI._define_camera_frame()

        # Trim the FOV and data array to the cutout window if necessary
        default_fov = inputs.fov
        if inputs.window is None:
            window = [0, 0, _SSI_PIXELS, _SSI_PIXELS]
        else:
            window = inputs.window
            if not full_fov:
                (l0, s0, l1, s1) = window
                default_fov = oops.fov.SliceFOV(default_fov, (s0, l0),
                                                (s1 - s0, l1 - l0))
                if data is not None:
                    data = data[l0:l1, s0:s1]

        # Define the Snapshot parameters
        params = {
            'cadence'     : oops.cadence.SnapCadence(tstart, inputs.texp),
            'fov'         : default_fov,
            'path'        : 'GLL',
            'frame'       : _SPICE_FRAME_NAME,
            'target'      : inputs.target,
            'lightsource' : inputs.lightsource,
            'calibrations': [],         # TBD
        }

        # Apply standard overrides
        oops.Host._apply_overrides(params, **overrides)

        obs = oops.observation.Snapshot(('v', 'u'), **params)
        if data is not None:
            obs.insert_subfield('data', data)

        # Extra attributes: use standardized names wherever practical
        obs.insert_subfield('dict', dict_)
        obs.insert_subfield('host', 'Galileo')
        obs.insert_subfield('instrument', 'SSI')
        obs.insert_subfield('time_from_sclk', inputs.tmid is None)
        for (name, key) in _ITEM_SUBFIELDS:
            if key in dict_:
                obs.insert_subfield(name, dict_[key])
        obs.insert_subfield('window', window)
        obs.insert_subfield('full_fov', full_fov)

        # With a custom frame, pointing never came from a CK, so any CK that happens to be
        # furnished is unrelated to this observation and must not be reported as used.
        kernels = _Galileo.used_kernels(obs.time, 'ssi',
                                        return_all_planets=return_all_planets,
                                        ck=using_cks)
        obs.insert_subfield('spice_kernels', kernels)
        obs.insert_subfield('spice_to_frame', oops.Matrix3.IDENTITY)
        obs.insert_subfield('spice_frame_name', _SPICE_FRAME_NAME)
        obs.insert_subfield('spice_frame_id', _SPICE_FRAME_ID)
        obs.insert_subfield('ck_object_id', _CK_OBJECT_ID)
        obs.insert_subfield('ck_sclk_id', _CK_SCLK_ID)
        obs.insert_subfield('spice_frame_is_frozen', False)

        return obs

    @staticmethod
    def _select(obs, select):
        """The choice that `select` makes from the one Snapshot of an image or a row.

        Parameters:
            obs (Snapshot): The observation.
            select (int | slice | str | tuple[int | slice | str, ...] | None): The choice,
                made by ``Host._apply_select`` of oops as from a tuple holding `obs`; None
                for `obs` itself.

        Returns:
            Snapshot | tuple[Snapshot, ...]: The choice.
        """

        if select is None:
            return obs
        return oops.Host._apply_select((obs,), select)

    @staticmethod
    def _pds3_options(pds3_options):
        """The PDS3 parser options of this host, updated by the caller's.

        Parameters:
            pds3_options (dict[str, Any] | None): The caller's options.

        Returns:
            dict[str, Any]: The options.
        """

        return {**_PDS3_DEFAULTS, **(pds3_options or {})}

    @staticmethod
    def tstart_from_sclk(label_dict):
        """The start of an SSI exposure, derived from its spacecraft clock count.

        The count names the frame in which the image was read out. The shutter opens a
        fixed delay after the start of that frame, which depends on FRAME_DURATION, or
        one frame earlier for an EXTENDED exposure taken before the phase-2 flight
        software of the Jupiter tour, whose first frame is GO_0017 C0346405900R (RIM
        3464059). The result agrees with the start of the exposure implied by IMAGE_TIME
        to within about 1/3 s for most images.

        The leap-second and spacecraft clock kernels are furnished on the first call, so
        no call to :meth:`initialize` is needed.

        Parameters:
            label_dict (collections.abc.Mapping[str, Any]): A PDS3 label of an SSI image,
                as a dictionary, or a row of an SSI index merged with its supplement. It
                must give SPACECRAFT_CLOCK_START_COUNT, in the label form "RRRRRRRR.MM"
                or the SPICE form "RRRRRRRR:MM:S:T", and FRAME_DURATION in seconds.
                EXPOSURE_TYPE, if absent, is taken as "NORMAL".

        Returns:
            float: The start of the exposure in seconds TDB.

        Raises:
            KeyError: If SPACECRAFT_CLOCK_START_COUNT or FRAME_DURATION is missing.
            oops.host.HostFormatError: If FRAME_DURATION is not one of the SSI frame
                durations, 2.333, 8.667, 15.167, 30.333, and 60.667 s, to within 0.05 s.
            oops.OopsValueError: If the clock count has more than four fields or a field
                that is not an unsigned integer.
        """

        (count, delay) = GalileoSSI._sclk_start(label_dict)
        return _Galileo.tdb_from_sclk(count) + delay

    @staticmethod
    def _sclk_start(dict_):
        """The clock count of an SSI image, and the delay of its shutter's opening.

        Parameters:
            dict_ (collections.abc.Mapping[str, Any]): The label or index row, as
                described for :meth:`tstart_from_sclk`.

        Returns:
            tuple[str, float]: The SPACECRAFT_CLOCK_START_COUNT, and the seconds from the
            start of the frame that it names to the opening of the shutter, negative for
            an EXTENDED exposure of phase 1.

        Raises:
            KeyError: If SPACECRAFT_CLOCK_START_COUNT or FRAME_DURATION is missing.
            oops.host.HostFormatError: If FRAME_DURATION is not one of the SSI frame
                durations.
        """

        count = dict_['SPACECRAFT_CLOCK_START_COUNT']
        frame_duration = float(dict_['FRAME_DURATION'])
        nominal = min(_SHUTTER_DELAY, key=lambda key: abs(key - frame_duration))
        if abs(nominal - frame_duration) > _FRAME_DURATION_TOLERANCE:
            raise oops.host.HostFormatError(f'unrecognized Galileo SSI FRAME_DURATION: '
                                            f'{frame_duration}')

        delay = _SHUTTER_DELAY[nominal]
        phase1 = int(_SCLK_SEPARATORS.split(count.strip())[0]) < _PHASE2_FIRST_RIM
        if dict_.get('EXPOSURE_TYPE', 'NORMAL') == 'EXTENDED' and phase1:
            delay -= nominal

        return (count, delay)

    ######################################################################################
    # Detectors
    ######################################################################################

    @staticmethod
    def _detect_in_pds3(label):
        """True if the given parsed PDS3 label describes this host's data."""

        # Labels spell the INSTRUMENT_NAME as "SOLID STATE IMAGING SYSTEM" or
        # "SOLID_STATE_IMAGING"
        return (label.get('SPACECRAFT_NAME', '').startswith('GALILEO')
                and label.get('INSTRUMENT_NAME', '').startswith('SOLID'))

    @staticmethod
    def _detect_in_vicar(label):
        """True if the given VicarLabel describes this host's data."""

        return label.get('MISSION', '') == 'GALILEO' and label.get('SENSOR', '') == 'SSI'

    @staticmethod
    def _detect_in_index(label):
        """True if the given index label refers to this host/instrument.

        The label of every current SSI metadata index names the host and instrument at
        its top level, so the rows are never needed.
        """

        return (label.get('INSTRUMENT_HOST_NAME', '').startswith('GALILEO')
                and label.get('INSTRUMENT_ID', '') == 'SSI')

    ######################################################################################
    # Initialization
    ######################################################################################

    @staticmethod
    def initialize(*, planets=None, asof=None, irregulars=True):
        """Initialize key information about the SSI instrument.

        Fills in key information about the camera and loads the SPICE kernels.
        :meth:`from_file` and :meth:`from_index` call it with the defaults, so it needs to
        be called first only to choose other options. After the first call, later calls
        to this function are ignored. Importing the module loads nothing.

        Parameters:
            planets (list, optional): A list of planets to pass to
                :meth:`~oops.Body.define_solar_system`. None or 0 means all.
            asof (str, optional): Only use SPICE kernels that existed before this date;
                None to ignore.
            irregulars (bool, optional): True to include the irregular satellites; False
                otherwise.

        Raises:
            oops.host.HostError: If the instrument kernel does not give an 800x800 image.
        """

        # Quick exit after first call
        if GalileoSSI._initialized:
            return

        _Galileo.initialize(planets=planets, asof=asof, irregulars=irregulars)
        _Galileo.load_instruments(asof=asof)

        # Load the instrument kernel
        GalileoSSI._INSTRUMENT_KERNEL = _Galileo.spice_instrument_kernel('SSI',
                                                                         asof=asof)[0]

        # Construct the FOVs
        info = GalileoSSI._INSTRUMENT_KERNEL['INS'][_NAIF_ID]

        cf = cspyce.gdpool(f'INS{_NAIF_ID}_DISTORTION_COEFF', 0)[0]
        focal = cspyce.gdpool(f'INS{_NAIF_ID}_FOCAL_LENGTH', 0)[0]
        px = cspyce.gdpool(f'INS{_NAIF_ID}_PIXEL_SIZE', 0)[0]
        cxy = cspyce.gdpool(f'INS{_NAIF_ID}_FOV_CENTER', 0)

        scale = px/focal
        distortion_coeff = [1, 0, cf / scale**2]

        shape = (info['MAX_SAMPLE'], info['MAX_LINE'])
        if shape != (_SSI_PIXELS, _SSI_PIXELS):
            raise oops.host.HostError(f'unsupported SSI image shape {shape} in the SSI '
                                      'kernel')

        fov_full = oops.fov.BarrelFOV(scale, shape,
                                      coefft_uv_from_xy=distortion_coeff,
                                      uv_los=(cxy[0], cxy[1]))
        fov_summed = oops.fov.SubsampledFOV(fov_full, 2)
#        fov_his =
#        fov_hma =
#        fov_hca = oops.fov.GapFOV(oops.fov.SubsampledFOV(fov_full, (1,4)),
#                                  (1,0.25))
#               ... maybe need SparseFOV or SkipFOV class
#        fov_him =

        # Construct FOV dictionary, keyed by the telemetry format ID
        fovs = {}
        fovs['FULL'] = fov_full

        # Phase-2 Telemetry Formats
        fovs['HIS'] = fov_summed
        fovs['HMA'] = fov_full
        fovs['HCA'] = fov_full
        fovs['HIM'] = fov_full
        fovs['IM8'] = fov_full
        fovs['AI8'] = fov_summed
        fovs['IM4'] = fov_full

        # Phase-1 Telemetry Formats
        fovs['XCM'] = fov_full
#        fovs['XED'] = fov_full
        fovs['HCJ'] = fov_full      # Inference based on inspection
        fovs['HCM'] = fov_full      # Inference based on inspection
                                    # hmmm, actually C0248807700R.img is 800x200
                                    # maybe this is just a cropped full fov
        fovs['NONE'] = fov_full     # Inference based on inspection

        GalileoSSI._FOVS = fovs

        # Load kernels
        _Galileo.load_kernels()

        # Define star HD 151288, a target not in the default catalog
        ra = 15*(16 + 45/60 + 6.3505410601/3600)    # 16 45 06.3505410601 in Simbad
        dec = 33 + 30/60 + 33.213816014/3600        # +33 30 33.213816014
        oops.lightsource.add_star('HD 151288', ra, dec, 8.11, 'K7.5Ve')

        GalileoSSI._initialized = True

    @staticmethod
    def _define_camera_frame():
        """Register the SPICE-derived GLL_SCAN_PLATFORM frame.

        :meth:`initialize` must have been called first. Built lazily (and only once) so
        that observations using a custom frame never construct or depend on the SPICE
        camera frame.
        """

        # The frame is tested against the Frame registry rather than against a flag on
        # this class. Registering a second frame under an ID that already exists does not
        # fail; it is renamed GLL_SCAN_PLATFORM_2 and never found again, so the test has
        # to ask the registry what it actually holds.
        if not oops.Frame.frame_id_exists(_SPICE_FRAME_NAME):
            _ = oops.frame.SpiceType1Frame(_SPICE_FRAME_NAME, _TICK_TOLERANCE)

    @staticmethod
    def _reset():
        """Reset the internal Galileo SSI parameters.

        Can be useful for debugging.
        """

        GalileoSSI._INSTRUMENT_KERNEL = None
        GalileoSSI._FOVS = {}
        GalileoSSI._initialized = False
        _Galileo.reset()


GalileoSSI._register()

##########################################################################################
