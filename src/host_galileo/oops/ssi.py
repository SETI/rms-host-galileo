##########################################################################################
# host_galileo/oops/ssi.py
##########################################################################################
"""The Galileo SSI host, which makes a Snapshot from each Solid State Imaging image.

The :doc:`user guide </user_guide/user_guide>` describes the products it reads, how it
merges the two labels of each product, and the errors it raises.
"""

import re

import julian
import cspyce

from filecache import FCPath
from pdsparser import Pds3Label
from vicar     import VicarError, VicarImage, VicarLabel

import oops
from ._galileo import _Galileo
from .         import _ssi_labels

__all__ = ['GalileoSSI']

# The SSI camera frame is the scan platform frame, whose C kernels are sampled at discrete
# spacecraft clock ticks. SSI images are spaced as closely as 1 unit in the file name,
# which corresponds to 80 clock ticks. Therefore, we use a tolerance of +/-40.
_SPICE_FRAME_NAME = 'GLL_SCAN_PLATFORM'
_SPICE_FRAME_ID = -77001
_TICK_TOLERANCE = 40

# The NAIF ID of the SSI instrument, as used in its instrument kernel
_NAIF_ID = -77036

# The telemetry format of an index row that gives none, as an index read without its
# supplement does; its FOV is the full frame
_INDEX_DEFAULT_MODE = 'NONE'

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
# frames of 30.333 s, but about 0.4 s earlier than normal at the other frame durations
# (0.721 s for 8.667 s, 0.775 s for 15.166 s, 6.068 s for 60.667 s). That difference is
# unused: all 77 frames with no IMAGE_TIME are phase-1 frames.
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

# Target name as found in PDS label -> correct name
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

# What the C kernels could not supply, as the OopsSpiceError for missing pointing names it
_NO_POINTING = 'no Galileo SSI pointing'


class GalileoSSI(oops.Host):
    """The Galileo SSI host: Snapshot constructors for Solid State Imaging images."""

    NAME = 'Galileo SSI'

    # Every SSI metadata index lists its files in the same order as its supplement, so the
    # rows are paired by position. A survey of all 23 of them (GO_0002-GO_0023 and the
    # cumulative GO_0999, 29,324 rows), reading VOLUME_ID and FILE_SPECIFICATION_NAME from
    # both tables, found the same files in the same order in every one, none listed twice.
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

        By default, the image and its field of view are trimmed to the label's cutout
        window.

        Given the VICAR image file, only that file is read. Given its PDS3 label, the
        image file that the label's IMAGE pointer gives, in the label's directory and in
        either upper or lower case, is read as well, unless `astrometry` is True. The
        image array is always located by the VICAR file's own header, never by the record
        number in the label's IMAGE pointer, so it is the same from either input.

        The metadata come from every label that was read: the PDS3 label and the VICAR
        header. A header value that the host converts into a PDS3 item, such as EXP into
        EXPOSURE_DURATION, is given only under its PDS3 name. A value that no label read
        gives is left out. The PDS3 null constants UNK, N/A, and NULL, and the header's
        unset values (a numeric -999, a PICNO or SPICE_C_ID of "?", a TCA of "CLOSEST
        TIME", a SMEAR of 0.1, and time fields that are not a valid time) count as no
        value, as do the texts UNKNOWN, NOT APPLICABLE, "---", and an empty string. So
        do the values that the archive's DATASET.CAT documents as not implemented, not
        applicable, or not available: a CONE_ANGLE (CONE), SMEAR_AZIMUTH (SMRAZ),
        CENTER_RING_RADIUS (RAD), CENTRAL_BODY_DISTANCE (PLRANGE), HORIZONTAL_PIXEL_SCALE
        (HSCL), VERTICAL_PIXEL_SCALE (VSCL), or SLANT_DISTANCE (SLRANGE) of 0; a
        SOLAR_DISTANCE (SOLRANGE) of 777909100 km; a time from closest approach that
        begins with "-73T"; and each of these groups when every value of it that a label
        gives is 0 and the label gives at least two: RIGHT_ASCENSION, DECLINATION, and
        TWIST_ANGLE (RA, DEC, and TWIST), as they are when the C kernels had no pointing;
        CENTER_LATITUDE and CENTER_LONGITUDE; INTERCEPT_POINT_LATITUDE and
        INTERCEPT_POINT_LONGITUDE (LAT and LON); and NORTH_AZIMUTH, SUB_SOLAR_AZIMUTH, and
        SUB_SPACECRAFT_AZIMUTH (NORAZ, SUNAZ, and SCAZ). One value of a group alone at 0
        is kept. A negative HORIZONTAL_PIXEL_SCALE (HSCL), VERTICAL_PIXEL_SCALE (VSCL), or
        SLANT_DISTANCE (SLRANGE), which is never real, counts as no value too. The time
        from closest approach, which the labels before the Jupiter
        tour call TIME_FROM_CLOSEST_APPROACH and SATELLITE_TIME_FROM_CLST_APR, is given
        under the tour labels' names, NTV_TIME_FROM_CLOSEST_APPROACH and
        NTV_SAT_TIME_FROM_CLOSEST_APRH, which the header's TCA also takes.

        Where both labels give a value and they differ, each item takes its value from
        one fixed label for every image: the VICAR header for EXPOSURE_DURATION and for
        the pointing and viewing geometry that the label rounds or miscomputes (such as
        DECLINATION, TARGET_CENTER_DISTANCE, and SUB_SOLAR_LONGITUDE) and for the
        compression ratios, and the PDS3 label for every other item, including
        TARGET_NAME. If the winning label lacks a value, the other label's value is
        used.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | Pds3Label | VicarLabel | VicarImage):
                The path to the PDS3 label or the VICAR image file, or its parsed content
                as a Pds3Label, VicarLabel, or VicarImage, which is read as that kind of
                input. A path ending in ".lbl", in any case, is a PDS3 label. The
                Snapshot of parsed content that has no file path lacks the path
                subfields; a full read of such a Pds3Label or VicarLabel, which needs the
                image file, raises an error. A VicarImage with a file path has its header
                read again from that file, so that `header_dict` holds the header as
                written; one without a file path, as VicarImage.from_file() gives, has
                only the header that it holds, whose HOST and BHOST rms-vicar sets to
                the reading machine.
            full_fov (bool, optional): If True, the full image is returned with a mask
                describing the regions with no data; otherwise, the image and FOV are
                trimmed to the label's cutout window.
            return_all_planets (bool, optional): Include kernels for all planets, not just
                the target of the mission phase.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser.
            astrometry (bool, optional):
                True to specify that the returned Observations will only be used for
                timing and astrometry. In this case, no data arrays are read and returned
                in the Observation.
            select (None): Must be None. An SSI image is a single observation, so there
                is nothing to select from; any other value raises an error.
            target (str | oops.Body, optional): Override for the default target Body of
                the observation. Use "NONE" for inertial pointing, indicating that no
                Solar System body was tracked.
            lightsource (oops.LightSource, optional): Override for the default LightSource
                of the observation.
            path (str | oops.Path, optional): Override for the Path of the observer.
            frame (str | oops.Frame, optional): Override for the Frame of the observing
                instrument.
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
            Snapshot: The observation, with `data` inserted unless `astrometry` is True,
            and these subfields:

            - `dict`: the merged metadata, keyed by the PDS3 names for the values that
              both labels give, plus every other value of each label read under its own
              name. A VICAR keyword the host converts appears only under its PDS3 name:
              RIM and MOD91 as SPACECRAFT_CLOCK_START_COUNT, the SCET fields as
              IMAGE_TIME. MOD10 and MOD8, the finer fields of the clock count, keep their
              own names. Times are ISO 8601 strings, such as "1996-06-26T08:46:04.890",
              sets are lists, and every value has the same type from either label, so
              that `json` can write the dictionary. The keys that the PDS3 parser adds
              (such as "^IMAGE_offset" and "IMAGE_TIME_fmt") are left out.
            - `label_dict`: the PDS3 label as given by `Pds3Label.as_dict()`, if the PDS3
              label was read, in a form that `json` can write: each time as the label
              writes it, each set as a list in label order, no key whose value is None,
              such as END, and none of the keys that the parser adds beside the
              keywords (such as "^IMAGE_offset" and "IMAGE_TIME_fmt").
            - `header_dict`: every parameter of the VICAR header as the file writes it,
              if it was read; a full read gives the same values as an astrometry read,
              including HOST and BHOST. A parameter that occurs more than once, such as
              the TASK of each step of the processing history, is keyed by its name and
              its occurrence counted from zero, as "TASK[0]", "TASK[1]", and so on.
            - `host` ("Galileo"), `instrument` ("SSI"), `time_from_sclk` (True if no
              IMAGE_TIME was read, so that the time comes from the spacecraft clock
              count), `window` (the cutout window as zero-based line and sample limits
              ``[l0, s0, l1, s1]``, or the full 800x800 frame if there is none), and
              `full_fov`.
            - `filter` (FILTER_NAME), `label_target` (TARGET_NAME), `observation_id`
              (OBSERVATION_ID), and `spacecraft_clock_start_count`
              (SPACECRAFT_CLOCK_START_COUNT, as "RRRRRRRR.MM"), each spelled as the PDS3
              label spells it, and each only if a label read gives it.
            - `spice_kernels`, `spice_to_frame`, `spice_frame_name`, and
              `spice_frame_id`.
            - `filepath`, `basename`, `abspath`, and `image_url`, which name the image
              file, also when `fileinfo` is its PDS3 label, and are absent when the
              parsed content given has no file path. `abspath` is given only when this
              call retrieved the image file to read its array, and is then the local path
              of the file or of its copy in the local cache. An astrometry read, which
              does not read the array, never gives it, nor does a VicarImage passed in,
              whose array was read before the call.

        Raises:
            oops.OopsSpiceError: If no `frame` is given and the C kernels lack the
                pointing of the scan platform at the time at which the frame is evaluated
                (the midtime of the exposure, or of the `cadence` override, or the start
                or end that a `tracker` selects, plus the seconds of a `timeshift` once if
                its letters include "T" and once more if they include "F"), or if the
                time must come from the spacecraft clock count and the clock kernel
                cannot convert it. The message of the first is
                "no Galileo SSI pointing at <time> for <image file>: <SPICE error code>",
                where <time> is that time without the seconds of an "F" timeshift, and
                the cause of either is the SPICE exception.
            TypeError: If a keyword argument is not one of those listed.
            oops.OopsValueError: If `select` is not None, or if `pds3_options` names a
                parser method other than fast, strict, and loose that the parser rejects,
                or another invalid option.
            oops.host.PdsHostError: If none of the fast, strict, and loose parsers can
                read the PDS3 label, or if it has no IMAGE pointer.
            oops.host.VicarHostError: If the VICAR header or image cannot be parsed,
                including an image file shorter than its header says, or if the image's
                array has more than one band.
            oops.host.HostFormatError: If no label read gives the exposure duration, the
                telemetry format, or a time (IMAGE_TIME, or the spacecraft clock count
                and frame duration), if the telemetry format is not one of SSI's, if the
                time must come from the spacecraft clock count and FRAME_DURATION is not
                one of the SSI frame durations, or if a full read is asked of a parsed
                Pds3Label or VicarLabel that has no file path, with a message that says
                no file path is known.
            oops.host.WrongHostError: If `fileinfo`, or the parsed label or image given,
                is not a Galileo SSI PDS3 label or VICAR image.
            FileNotFoundError: If the image file that the PDS3 label names does not
                exist, unless `astrometry` is True.
            OSError: If a file cannot be read.
        """

        oops.Host._check_select(select, unit='Galileo SSI image')

        GalileoSSI.initialize()  # define everything the first time through; use defaults
                                 # unless initialize() was called explicitly.

        # Default Pds3Label options are method="fast"; the parser's source text is never
        # wanted, in label_dict or anywhere else
        pds3_options = oops.Host._pds3_options_with_defaults(pds3_options)
        pds3_options['sources'] = False

        fileinfo = GalileoSSI._parse_fileinfo(fileinfo, pds3_options=pds3_options)

        # Read the image file unless this is an astrometry read; never read the label
        # beside an image file. The header comes from a VicarLabel, never from a
        # VicarImage, which replaces the HOST and BHOST of the file with the reading
        # machine's; a VicarImage passed in with a file path has its header read again
        # from the file. This stays until rms-vicar keeps the header as the file writes
        # it (SETI/rms-vicar#34). Only a read that retrieves the image file here, for its
        # array, knows a local path that exists, so only it gives abspath.
        label_dict = None
        header = None
        data = None
        abspath = None
        if isinstance(fileinfo, Pds3Label):
            label_dict = _ssi_labels.label_json(fileinfo.as_dict())
            filepath = GalileoSSI._image_path(fileinfo, find=not astrometry)
            if not astrometry:
                filepath = GalileoSSI._known_path(filepath, 'PDS3 label')
                header = GalileoSSI._read_vicar_header(filepath)
                (data, abspath) = GalileoSSI._read_vicar_image(filepath)
        else:
            filepath = None if fileinfo.filepath is None else FCPath(fileinfo.filepath)
            if isinstance(fileinfo, VicarImage):
                header = (fileinfo.label if filepath is None
                          else GalileoSSI._read_vicar_header(filepath))
                if not astrometry:
                    data = GalileoSSI._array_2d(fileinfo, filepath)
            else:
                header = fileinfo
                if not astrometry:
                    filepath = GalileoSSI._known_path(filepath, 'VICAR header')
                    (data, abspath) = GalileoSSI._read_vicar_image(filepath)

        header_dict = None if header is None else _ssi_labels.header_json(header)

        dict_ = _ssi_labels.merge_items(
            None if label_dict is None else _ssi_labels.label_items(label_dict),
            None if header_dict is None else _ssi_labels.header_items(header_dict))

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

        snapshot = GalileoSSI._make_snapshot(dict_, filepath=filepath, abspath=abspath,
                                             data=data, full_fov=full_fov,
                                             return_all_planets=return_all_planets,
                                             overrides=overrides, check_pointing=True)
        if label_dict is not None:
            snapshot.insert_subfield('label_dict', label_dict)
        if header_dict is not None:
            snapshot.insert_subfield('header_dict', header_dict)
        return snapshot

    @staticmethod
    def _parse_fileinfo(fileinfo, *, pds3_options):
        """The parsed PDS3 label or VICAR header of a Galileo SSI image.

        A PDS3 label that the requested parser rejects, whether for its syntax or for a
        value it cannot convert (such as a time with an hour of 25), is tried again with
        each later parser of fast, strict, and loose, as oops.Host._read_fileinfo() does
        with `pds3_fallback`. A VICAR file's header is read without its array.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | Pds3Label | VicarLabel | VicarImage):
                The path to the PDS3 label or the VICAR image file, or its parsed content.
            pds3_options (dict[str, Any]): The options for the PDS3 label parser,
                including its method.

        Returns:
            Pds3Label | VicarLabel | VicarImage: The parsed content.

        Raises:
            oops.OopsValueError: If the requested method is none of fast, strict, and
                loose and the parser rejects it or another option.
            oops.host.PdsHostError: If no parser accepts the PDS3 label. The message
                names the label and gives the error of each parser tried.
            oops.host.VicarHostError: If the VICAR header cannot be parsed.
            oops.host.WrongHostError: If the file or the parsed content is not that of a
                Galileo SSI image.
            OSError: If the file cannot be read.
        """

        parsed = GalileoSSI._read_fileinfo(fileinfo, formats='PV', astrometry=True,
                                           pds3_options=pds3_options, pds3_fallback=True)

        # A parsed label or image passed in has not been identified yet
        if isinstance(parsed, Pds3Label):
            if not GalileoSSI._detect_in_pds3(parsed):
                raise oops.Host._no_host_error([GalileoSSI], 'PDS3 label',
                                               parsed.filepath)
        elif not GalileoSSI._detect_in_vicar(parsed):
            raise oops.Host._no_host_error([GalileoSSI], 'VICAR file', parsed.filepath)

        return parsed

    @staticmethod
    def _image_path(label, *, find):
        """The image file that a PDS3 label describes.

        Parameters:
            label (Pds3Label): The parsed label.
            find (bool): True to look for the file, which then need not have the same
                case as in the label; False to name it without touching the file system.

        Returns:
            FCPath | None: The file named by the label's IMAGE pointer, in the directory
            of the label. If `find` is True and no file of that name exists, the same
            name in upper or lower case is used if it exists. None if the label has no
            file path.

        Raises:
            oops.host.PdsHostError: If the label has no IMAGE pointer.
        """

        name = label.get('^IMAGE')
        if not isinstance(name, str):
            raise oops.host.PdsHostError(f'no ^IMAGE pointer in Galileo SSI PDS3 label '
                                         f'{label.filepath}')
        if label.filepath is None:
            return None
        parent = FCPath(label.filepath).parent
        if find:
            for option in dict.fromkeys((name, name.upper(), name.lower())):
                if (parent / option).is_file():
                    return parent / option
        return parent / name

    @staticmethod
    def _known_path(filepath, kind):
        """The path to the image file, which a full read needs.

        Parameters:
            filepath (FCPath | None): The image file, or None if no path is known.
            kind (str): What was given, such as "PDS3 label", named in the error.

        Returns:
            FCPath: `filepath`.

        Raises:
            oops.host.HostFormatError: If `filepath` is None.
        """

        if filepath is None:
            raise oops.host.HostFormatError(f'no file path is known for the image file '
                                            f'of this Galileo SSI {kind}')
        return filepath

    @staticmethod
    def _read_vicar_header(filepath):
        """The VICAR header of an image file, as the file writes it.

        Parameters:
            filepath (FCPath): The image file.

        Returns:
            VicarLabel: The header.

        Raises:
            oops.host.VicarHostError: If the file has no valid VICAR header.
            OSError: If the file cannot be read.
        """

        try:
            header = VicarLabel(filepath, strict=False)
        except VicarError as err:
            raise oops.host.VicarHostError(f'{err}: {filepath}') from err
        return header

    @staticmethod
    def _read_vicar_image(filepath):
        """The full array of a VICAR image file, and the local path it was read from.

        The file is retrieved first, so that a remote file has a local copy whose path
        exists when it is returned.

        Parameters:
            filepath (FCPath): The image file.

        Returns:
            tuple[numpy.ndarray, pathlib.Path]: The array, located by the file's own
            header, and the local path that ``FCPath.retrieve()`` returned.

        Raises:
            oops.host.VicarHostError: If the file is not a valid VICAR image, including
                one whose array is shorter than its header says, or if its array has
                more than one band.
            OSError: If the file cannot be retrieved or read.
        """

        local = FCPath(filepath).retrieve()

        # A truncated array fails to take the header's shape with a bare ValueError
        try:
            vicar = VicarImage(local, strict=False)
        except (VicarError, ValueError) as err:
            raise oops.host.VicarHostError(f'{err}: {filepath}') from err
        return (GalileoSSI._array_2d(vicar, filepath), local)

    @staticmethod
    def _array_2d(image, filepath):
        """The array of a VICAR image, in two dimensions.

        Parameters:
            image (VicarImage): The image.
            filepath (FCPath | None): The image file, for the message; None if the image
                has no file path.

        Returns:
            numpy.ndarray: The array.

        Raises:
            oops.host.VicarHostError: If the array has more than one band.
        """

        try:
            return image.data_2d
        except VicarError as err:
            raise oops.host.VicarHostError(str(err) if filepath is None
                                           else f'{err}: {filepath}') from err

    @staticmethod
    def from_index(filepath, *, supplement=oops.Host.SUPPLEMENT, label=None,
                   row_dicts=None, pds3_options=None, filter=None, full_fov=False,
                   return_all_planets=False, select=None, target=None, lightsource=None,
                   path=None, frame=None, cadence=None, fov=None, calibrations=None,
                   timeshift=None, navigation=None, tracker=None, parallel=None):
        """An iterator over Snapshots, one for each row of a Galileo SSI index file.

        It reads a Ring-Moon Systems Node metadata index and its supplement. The call
        itself checks the arguments and reads the labels of the index and its
        supplement, so that a wrong call fails at once. The rows of the two tables are
        then read as the iteration proceeds, a few at a time, and each row's Snapshot is
        built when the iteration reaches that row. The memory used stays the same
        however long the index is only with an rms-oops that frees the frames and paths
        built for each row once its Snapshot is gone; with one that keeps them, it grows
        a little with each row. An error in a row is raised when the iteration
        reaches it, after every earlier row has been yielded (except for a row longer
        than its label says, as described below), and the iteration ends there;
        `filter` can leave such a row out.

        Every SSI index lists its files in the same order as its supplement, so the two
        tables are read in step and each index row is paired with the supplement row in
        the same position. The two rows must give the same FILE_SPECIFICATION_NAME.

        Parameters:
            filepath (str | pathlib.Path | FCPath): The full path to the PDS label of a
                Galileo SSI index of the Ring-Moon Systems Node metadata tree, such as
                ".../metadata/GO_0xxx/GO_0017/GO_0017_index.lbl".
            supplement (str | pathlib.Path | FCPath | object, optional): The path to the
                PDS label of the supplemental index file, whose columns are merged into
                each row; an SSI index needs its supplement. Where both tables give a
                column, the index value is kept. The default, `oops.Host.SUPPLEMENT`,
                reads the file in the directory of `filepath` whose name matches except
                that "_index" is replaced by "_supplemental_index" in upper or lower
                case.
            label (Pds3Label, optional): The PDS label of `filepath` if it has already
                been read.
            row_dicts (collections.abc.Iterable[dict[str, Any]], optional): The rows of
                the index, each already merged with its supplement row, if the tables
                have already been read. They are iterated over in place of the tables,
                and `supplement` is not used. The label of `filepath` is still read,
                unless `label` is given, because it tells whether the table writes the
                compression ratios in rotated columns.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser.
            filter (dict[str, Any], optional): A dictionary of index column names and
                their values or a set or tuple of values. If provided, only rows of the
                index in which each named column equals the value, falls in the set of
                values, or falls in a range defined by a tuple of two values, will be
                included. A column that neither table has raises an error, and a row
                without a named column is left out.
            full_fov (bool, optional): If True, each field of view covers the full image
                rather than the cutout window.
            return_all_planets (bool, optional): Include kernels for all planets, not just
                the target of the mission phase.
            select (None): Must be None. Each row describes a single observation, so
                there is nothing to select from; any other value raises an error.
            target (str | oops.Body, optional): Override for the default target Body of
                every observation. Use "NONE" for inertial pointing, indicating that no
                Solar System body was tracked.
            lightsource (oops.LightSource, optional): Override for the default LightSource
                of every observation.
            path (str | oops.Path, optional): Override for the Path of the observer.
            frame (str | oops.Frame, optional): Override for the Frame of the observing
                instrument.
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
            collections.abc.Iterator[Snapshot]: One observation for every row of the
            index that passes `filter`, in index order, including any row with a zero
            exposure, each with subfields `dict`, `host`, `instrument`, `filter`,
            `time_from_sclk`, `label_target`, `observation_id`, `spice_kernels`,
            `filepath`, `basename`, `image_url`, `spice_to_frame`, `spice_frame_name` and
            `spice_frame_id` inserted, and `spacecraft_clock_start_count` if the row gives
            it. No row has `abspath`, because the image files are not retrieved.
            `dict` is the row, merged with its row of the supplement, with every value as
            a Python str, int, float, or list. A value that stands for "no value" is left
            out: a PDS3 null such as UNK, a numeric -999, the placeholder a column writes
            for an image that has no such value (the supplement's CUT_OUT_WINDOW and
            TRUTH_WINDOW of -1, -1, -1, -1, its ICT_DESPIKE_THRESHOLD and
            ICT_QUANTIZATION_STEP_SIZE of -1, a compression ratio of 0, and a
            SMEAR_MAGNITUDE of 0.1), the values that the archive documents as standing
            for no value, by the same rules as for a label read by :meth:`from_file`
            (such as a SOLAR_DISTANCE of 777909100 km, a HORIZONTAL_PIXEL_SCALE or
            VERTICAL_PIXEL_SCALE of 0 or less, a negative SLANT_DISTANCE, a time from
            closest approach that begins with
            "-73T", and a RIGHT_ASCENSION, DECLINATION, and TWIST_ANGLE that are all 0),
            and the ORBIT_NUMBER of -99 that the index gives a frame taken before the
            Jupiter tour. Three values that the tables declare as constants are kept,
            because they are real: an EXPOSURE_DURATION of 0, an ON_CHIP_MOSAIC_FLAG of
            N, and a PRODUCT_VERSION_ID of 0. Where the index label names the three
            compression-ratio columns in the order of their bytes that every SSI index
            label has used, the values are given under the names of the ratios they hold,
            as the product labels give them. The DATA_SET_ID of an
            image of the Jupiter tour has the version 1.1 of the product labels, where
            the current indexes write 1.0. NTV_SAT_TIME_FROM_CLOSEST_APR is given under
            the labels' name, NTV_SAT_TIME_FROM_CLOSEST_APRH.
            The C kernels are checked for pointing when the Snapshot is built, but a
            failed check does not raise: for an image whose pointing they cannot supply,
            the default frame raises oops.OopsSpiceError when it is evaluated. A `frame`
            override is not checked.
            `filepath`, `basename`, and `image_url` name the image file, as a
            full path or URL: the row's FILE_SPECIFICATION_NAME, with the label's
            extension ".LBL" replaced by ".IMG" in the same case, in the volume directory
            named by its VOLUME_ID. The volume directories are found from the location of
            `filepath`: under "volumes" in place of "metadata" for an index in an RMS
            Node metadata tree (".../metadata/GO_0xxx/GO_0017/GO_0017_index.lbl"), and
            otherwise in the index's own directory.

        Raises:
            TypeError: At call time, if a keyword argument is not one of those listed.
            oops.OopsTypeError: When the iteration reaches a row, if a range in `filter`
                cannot be compared with the row's value.
            oops.OopsValueError: At call time, if `select` is not None, if a value in
                `pds3_options` is invalid, or if `filter` uses a column that neither
                table has. With `row_dicts`, `filter` is checked against the first row:
                at call time if `row_dicts` is a sequence, and when the iteration reads
                the first row if it is an iterator.
            oops.host.HostIndexError: At call time, if `supplement` is None or empty
                while `row_dicts` is not given; if `supplement` is
                `oops.Host.SUPPLEMENT` and the name of `filepath` contains no "_index";
                if either label does not describe a table that can be read; or if a
                table lacks FILE_SPECIFICATION_NAME. When the iteration
                reaches a row, if the index row and the supplement row in the same
                position give different files, or if the supplement ends before the
                index; at the end of a table file that holds fewer rows than its label
                says, after every row it holds; and at the first of the 100 rows read
                with a row longer than the label's row length, none of which is yielded.
                After the last index row, if the supplement has rows left over or a
                table file holds more rows than its label says. A value that cannot be
                read as its column's number is given as a string, with a UserWarning
                from the table reader, and the Snapshot of its row may then raise the
                TypeError or ValueError of the value's use.
            oops.host.HostFormatError: When the iteration reaches a row that lacks the
                exposure duration or a time (IMAGE_TIME, or the spacecraft clock count
                and frame duration), or whose time must come from the spacecraft clock
                count while its FRAME_DURATION is not one of the SSI frame durations.
            oops.OopsSpiceError: When the iteration reaches a row whose time must come
                from the spacecraft clock count and the clock kernel cannot convert it.
            oops.host.PdsHostError: At call time, if the label of the index or of the
                supplement cannot be parsed.
            oops.host.WrongHostError: At call time, if `filepath` is not a Galileo SSI
                index.
            OSError: At call time, if the label of `filepath` or of the supplement does
                not exist or cannot be read (FileNotFoundError for a missing
                supplement). When the iteration starts, if a table file does not exist or
                cannot be read; at call time for the index table, if the index label
                lacks INSTRUMENT_HOST_NAME and INSTRUMENT_ID, so that the first row is
                read to recognize the index.
        """

        oops.Host._check_select(select, unit='row of a Galileo SSI index')

        # The messages name the index as an FCPath, as oops.Host.from_index passes it, so
        # that they read the same on every platform and from either entry point
        index_path = FCPath(filepath)

        # The supplement is required, so None and the empty path are refused rather than
        # taken to mean "read the index alone", as oops takes them
        if row_dicts is None and (supplement is None or supplement == ''):
            raise oops.host.HostIndexError(f'Galileo SSI index needs its supplemental '
                                           f'index: {index_path}')

        GalileoSSI.initialize()

        # The index label decides, once for the whole table, whether the compression
        # ratios are in rotated columns, so it is read also when the rows are given. It
        # is parsed once and passed on.
        label = GalileoSSI._read_index_label(index_path, label=label,
                                             pds3_options=pds3_options)
        rotated = _ssi_labels.ratios_rotated(label)

        # Read the index and the supplement as the iteration proceeds, a few rows at a
        # time, each index row paired with its supplement row and filtered. The call
        # itself reads only the labels. Both tables also give
        # SPACECRAFT_CLOCK_START_COUNT, but the supplement writes it in the SPICE form
        # "RRRRRRRR:MM:S:T" and the index in the label form "RRRRRRRR.MM", so the default
        # validation, which compares every shared column, would reject every real pair
        # until the index is corrected (rms-data-projects #258, item 4); the pairing's
        # comparison of FILE_SPECIFICATION_NAME is the only check.
        rows = GalileoSSI._index_rows(index_path, supplement=supplement, label=label,
                                      row_dicts=row_dicts, pds3_options=pds3_options,
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

        # Build each row's Snapshot when the iteration reaches it.
        # FILE_SPECIFICATION_NAME names the label; the Snapshot describes the image beside
        # it.
        volumes = oops.Host._volumes_directory(index_path)

        def snapshots():
            for row_dict in rows:
                label_path = oops.Host._index_file_path(volumes, row_dict)
                (stem, _, suffix) = label_path.name.rpartition('.')
                image = label_path.with_name(stem + '.'
                                             + ('img' if suffix.islower() else 'IMG'))
                yield GalileoSSI._make_snapshot(
                    _ssi_labels.row_items(row_dict, rotated_ratios=rotated),
                    filepath=image, data=None, full_fov=full_fov,
                    return_all_planets=return_all_planets, overrides=overrides,
                    default_mode=_INDEX_DEFAULT_MODE)

        return snapshots()

    @staticmethod
    def _make_snapshot(dict_, filepath, data=None, *, abspath=None, full_fov=False,
                       return_all_planets=False, overrides, default_mode=None,
                       check_pointing=False):
        """A Snapshot built from the merged metadata of an image or a row of an index.

        Parameters:
            dict_ (dict): The merged metadata of the labels read by :meth:`from_file`, or
                the items of an index row from `_ssi_labels.row_items()`; either holds no
                null value.
            filepath (FCPath | None): The path to the image file, as a full path or URL
                for an index row; None if no path is known, to leave out the path
                subfields.
            data (numpy.ndarray, optional): The full image array; None to omit the `data`
                subfield.
            abspath (pathlib.Path, optional): The local path of the image file, given only
                when this call retrieved the file to read its array; None to leave out
                the `abspath` subfield.
            full_fov (bool, optional): True to keep the full FOV and image; False to trim
                both to the cutout window.
            return_all_planets (bool, optional): Include kernels for all planets, not just
                the target of the mission phase.
            overrides (dict[str, Any]): The values of the standard overrides, `target`,
                `lightsource`, `path`, `frame`, `cadence`, `fov`, `calibrations`,
                `timeshift`, `navigation`, `tracker` and `parallel`, as described for
                :meth:`from_file`.
            default_mode (str, optional): The telemetry format to use if `dict_` gives
                none; None to require one.
            check_pointing (bool, optional): True to require that the C kernels give the
                scan platform's pointing at the time at which the frame is evaluated, as
                described for :meth:`from_file`, unless a frame override is given; False
                to give the Snapshot a frame that raises oops.OopsSpiceError when it is
                evaluated if they do not.

        Returns:
            Snapshot: The observation, with the subfields listed for :meth:`from_file`.

        Raises:
            oops.host.HostFormatError: If `dict_` lacks the exposure duration, the
                telemetry format, or a time (IMAGE_TIME, or the spacecraft clock count
                and frame duration), if the telemetry format is not one of SSI's, or if
                the time must come from the spacecraft clock count and FRAME_DURATION is
                not one of the SSI frame durations.
            oops.OopsSpiceError: If `check_pointing` is True, no frame override is given,
                and the C kernels lack the pointing of the scan platform at the time
                described for :meth:`from_file`; or if the time must come from the
                spacecraft clock count and the clock kernel cannot convert it.
        """

        # PDS3 labels and index rows use the same names
        value = dict_.get

        def required(key, needed_for):
            item = value(key)
            if item is None:
                raise oops.host.HostFormatError(f'Galileo SSI {needed_for} unknown: no '
                                                f'{key} for {filepath}')
            return item

        texp = required('EXPOSURE_DURATION', 'exposure') / 1000.

        # IMAGE_TIME is the middle of the exposure. Some RAW_CAL frames in GO_0002 and
        # GO_0003 have no IMAGE_TIME, but every label carries the spacecraft clock count.
        # Never fall back to a placeholder time: it silently places the frame at the
        # J2000 epoch.
        image_time = value('IMAGE_TIME')
        if image_time is None:
            tstart = GalileoSSI._tstart_from_sclk(
                required('SPACECRAFT_CLOCK_START_COUNT', 'time'),
                required('FRAME_DURATION', 'time'), value('EXPOSURE_TYPE'), filepath)
            time_from_sclk = True
        else:
            tmid = julian.tdb_from_tai(julian.tai_from_iso(image_time))
            tstart = tmid - texp / 2.
            time_from_sclk = False

        mode = value('TELEMETRY_FORMAT_ID') or default_mode
        if mode is None:
            required('TELEMETRY_FORMAT_ID', 'telemetry format')
        if mode not in GalileoSSI._FOVS:
            raise oops.host.HostFormatError(f'unrecognized Galileo SSI '
                                            f'TELEMETRY_FORMAT_ID {mode!r} '
                                            f'for {filepath}')

        label_target = value('TARGET_NAME')
        target = _TARGET_NAME_REPAIRS.get(label_target, label_target)
        lightsource = 'SUN'
        if target is not None and _TARGET_STARS.match(target):
            lightsource = oops.lightsource.star_lookup(target)
            target = 'NONE'

        # Make sure the SPICE frame is defined unless a custom frame was given
        using_cks = overrides['frame'] is None
        if using_cks:
            GalileoSSI._define_camera_frame()

        # The cutout window is (line, sample, lines, samples), with a one-based origin
        window = value('CUT_OUT_WINDOW')

        # Define the cutout window; Trim the FOV and data array if necessary
        default_fov = GalileoSSI._FOVS[mode]
        if window is None:
            window = [0, 0, _SSI_PIXELS, _SSI_PIXELS]
        else:
            l1, s1, nl, ns = window
            window = [l1 - 1, s1 - 1, l1 + nl - 1, s1 + ns - 1]
            if not full_fov:
                l0, s0, l1, s1 = window
                default_fov = oops.fov.SliceFOV(default_fov, (s0, l0), (ns, nl))
                if data is not None:
                    data = data[l0:l1, s0:s1]

        # Define the Snapshot parameters
        params = {
            'cadence'     : oops.cadence.SnapCadence(tstart, texp),
            'fov'         : default_fov,
            'path'        : 'GLL',
            'frame'       : _SPICE_FRAME_NAME,
            'target'      : target,
            'lightsource' : lightsource,
            'calibrations': [],         # TBD
        }

        # Check the pointing at the time at which the Snapshot will evaluate its frame,
        # which a cadence override and a timeshift of the cadence ("T") or of the frame
        # ("F") can change. An index row without its pointing still gives a Snapshot,
        # whose frame raises the error when it is used.
        if using_cks:
            params['frame'] = oops.Host._check_spice_frame(
                params['frame'], params['cadence'], overrides,
                description=_NO_POINTING, filepath=filepath, defer=not check_pointing)

        # Apply standard overrides, including the path subfields
        oops.Host._apply_overrides(params, filepath, abspath=abspath, **overrides)

        obs = oops.observation.Snapshot(('v', 'u'), **params)
        if data is not None:
            obs.insert_subfield('data', data)

        # Extra attributes: use standardized names wherever practical
        obs.insert_subfield('dict', dict_)
        obs.insert_subfield('host', 'Galileo')
        obs.insert_subfield('instrument', 'SSI')
        obs.insert_subfield('time_from_sclk', time_from_sclk)
        for (name, key) in _ITEM_SUBFIELDS:
            if value(key) is not None:
                obs.insert_subfield(name, value(key))
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

        return obs

    @staticmethod
    def _tstart_from_sclk(sclk, frame_duration, exposure_type, filepath):
        """The start of the exposure derived from the spacecraft clock count.

        The count names the frame in which the image was read out. The shutter opens a
        fixed delay after the start of that frame, which depends on FRAME_DURATION, or
        one frame earlier for an EXTENDED exposure taken before the Jupiter tour's
        phase-2 flight software. The result agrees with the start of the exposure implied
        by IMAGE_TIME to within about 1/3 s for most images.

        Parameters:
            sclk (str): The SPACECRAFT_CLOCK_START_COUNT, as "RRRRRRRR.MM" or in the
                SPICE form "RRRRRRRR:MM:S:T".
            frame_duration (float | str): The FRAME_DURATION in seconds.
            exposure_type (str | None): The EXPOSURE_TYPE; None is taken as "NORMAL".
            filepath (FCPath | None): The image file, named in an error; None if it is
                not known.

        Returns:
            float: The start of the exposure in seconds TDB.

        Raises:
            oops.host.HostFormatError: If FRAME_DURATION is not one of the SSI frame
                durations.
            oops.OopsSpiceError: If the clock kernel cannot convert the count.
        """

        frame_duration = float(frame_duration)
        nominal = min(_SHUTTER_DELAY, key=lambda key: abs(key - frame_duration))
        if abs(nominal - frame_duration) > _FRAME_DURATION_TOLERANCE:
            raise oops.host.HostFormatError(f'unrecognized Galileo SSI FRAME_DURATION '
                                            f'{frame_duration} for {filepath}')

        tstart = _Galileo.tdb_from_sclk(sclk, filepath=filepath) + _SHUTTER_DELAY[nominal]
        phase1 = int(_SCLK_SEPARATORS.split(sclk.strip())[0]) < _PHASE2_FIRST_RIM
        if exposure_type == 'EXTENDED' and phase1:
            tstart -= nominal

        return tstart

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

        The label of an RMS Node metadata index names the host and instrument at its top
        level. For a label without the host's name this detector returns None, and the
        host is judged from a row by :meth:`_detect_in_row`.
        """

        host = label.get('INSTRUMENT_HOST_NAME')
        if host is not None:
            return host.startswith('GALILEO') and label.get('INSTRUMENT_ID', '') == 'SSI'
        return None

    @staticmethod
    def _detect_in_row(row_dict):
        """True if the given row of an index file label refers to this host/instrument."""

        return (row_dict.get('MISSION_NAME', '') == 'GALILEO'
                and row_dict.get('INSTRUMENT_ID', '') == 'SSI')

    ######################################################################################
    # Initialization
    ######################################################################################

    @staticmethod
    def initialize(*, planets=None, asof=None, irregulars=True):
        """Initialize key information about the SSI instrument.

        Fills in key information about the camera and loads the SPICE kernels.
        :meth:`from_file` and :meth:`from_index` call it with the defaults, so it needs to
        be called first only to choose other options. After the first call, later calls
        to this function are ignored.

        Parameters:
            planets (int | collections.abc.Sequence[int], optional): The number of one
                planet, from 1 for Mercury to 9 for Pluto, or a sequence of such numbers,
                such as (5,) for Jupiter alone; each planet and its moons are defined
                through :meth:`~oops.Body.define_solar_system`. None or 0 means all.
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
