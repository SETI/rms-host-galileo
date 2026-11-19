##########################################################################################
# host_galileo/oops/ssi.py
##########################################################################################

import re

import julian
import cspyce

from filecache import FCPath
from pdsparser import Pds3Label
from vicar     import VicarError, VicarImage

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

# The columns that identify the same image in an index and its supplement. Both tables
# also give SPACECRAFT_CLOCK_START_COUNT, but the supplement writes it in the SPICE form
# "RRRRRRRR:MM:S:T" and the index in the label form "RRRRRRRR.MM", so the default
# validation of oops.Host._read_index_rows(), which compares every shared column, rejects
# every real pair.
_INDEX_KEYS = ['VOLUME_ID', 'FILE_SPECIFICATION_NAME']

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

# The subfields that give a metadata item unchanged, and the items' names. Each subfield
# is inserted only if the item has a value.
_ITEM_SUBFIELDS = [
    ('filter'                      , 'FILTER_NAME'                 ),
    ('label_target'                , 'TARGET_NAME'                 ),
    ('observation_id'              , 'OBSERVATION_ID'              ),
    ('spacecraft_clock_start_count', 'SPACECRAFT_CLOCK_START_COUNT'),
]

# The PDS3 label parsers to try, in order. The fast parser rejects some valid labels, such
# as one whose quoted value starts on the line after its keyword.
_PDS3_METHODS = ('fast', 'strict', 'loose')

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

# The exceptions that cspyce raises for SPICE data the loaded kernels do not hold
_SPICE_ERRORS = (LookupError, OSError, RuntimeError, ValueError)


class GalileoSSI(oops.Host):
    """The Galileo SSI host: Snapshot constructors for Solid State Imaging images."""

    NAME = 'Galileo SSI'

    _INSTRUMENT_KERNEL = None
    _FOVS = {}
    _initialized = False

    @staticmethod
    def from_file(fileinfo, *, full_fov=False, return_all_planets=False,
                  pds3_options=None, astrometry=False, target=None, lightsource=None,
                  path=None, frame=None, cadence=None, fov=None, calibrations=None,
                  timeshift=None, navigation=None, tracker=None, parallel=None, **kwargs):
        """A Snapshot object based on a given Galileo SSI image file.

        By default, only the valid image region is returned.

        Given the VICAR image file, only that file is read. Given its PDS3 label, the
        image file that the label's IMAGE pointer names, in the label's directory and in
        either upper or lower case, is read as well, unless `astrometry` is True. The
        image array is always located by the VICAR file's own header, never by the record
        number in the label's IMAGE pointer, so it is the same from either input.

        The metadata come from every label that was read: the PDS3 label and the VICAR
        header, whose values the host converts into PDS3 items (such as EXP into
        EXPOSURE_DURATION) are given only under their PDS3 names. A value that no label
        read gives is left out; the PDS3 null constants UNK, N/A, and NULL, and the
        header's unset values (a numeric -999, a PICNO or SPICE_C_ID of "?", time fields
        that make no valid time), count as no value. Where both labels give a value and
        they differ, each item takes its value from one fixed label for every image: the
        VICAR header for EXPOSURE_DURATION and for the pointing and viewing geometry that
        the label rounds or miscomputes (such as DECLINATION, TARGET_CENTER_DISTANCE, and
        SUB_SOLAR_LONGITUDE) and for the compression ratios, and the PDS3 label for every
        other item, including TARGET_NAME. If the winning label gives no value, the other
        label's value is used.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | Pds3Label | VicarLabel | VicarImage):
                The path to the PDS3 label or the VICAR image file, or its parsed content
                as a Pds3Label, VicarLabel, or VicarImage, which is read as that kind of
                input. A path ending in ".lbl", in any case, is a PDS3 label.
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
            **kwargs: Additional keyword arguments; they are accepted and ignored. This
                includes `select`.

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
              writes it, each set as a list in label order, and no key whose value is
              None, such as END.
            - `header_dict`: every parameter of the VICAR header, if it was read. A
              parameter that occurs more than once, such as the TASK of each step of the
              processing history, is keyed by its name and its occurrence counted from
              zero, as "TASK[0]", "TASK[1]", and so on.
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
              file, also when `fileinfo` is its PDS3 label. `abspath` is absent after an
              astrometry read of a remote PDS3 label, which never fetches the image,
              unless an earlier read put the image in the local cache.

        Raises:
            oops.OopsSpiceError: If no `frame` is given and the C kernels hold no
                pointing for the scan platform at the midtime of the exposure (after any
                `cadence` override). Its message is "no Galileo SSI pointing at <time> for
                <image file>: <SPICE error code>", and its cause is the SPICE exception.
            oops.OopsValueError: If `pds3_options` names a parser method other than fast,
                strict, and loose that the parser rejects, or another invalid option.
            oops.host.PdsHostError: If none of the fast, strict, and loose parsers can
                read the PDS3 label, or if it has no IMAGE pointer.
            oops.host.VicarHostError: If the VICAR header or image cannot be parsed,
                including an image file shorter than its header says.
            oops.host.HostFormatError: If no label read gives the exposure duration, the
                telemetry format, or a time (IMAGE_TIME, or the spacecraft clock count
                and frame duration), if the telemetry format is not one of SSI's, or if
                the time must come from the spacecraft clock count and FRAME_DURATION is
                not one of the SSI frame durations.
            oops.host.WrongHostError: If `fileinfo`, or the parsed label or image given,
                is not a Galileo SSI PDS3 label or VICAR image.
            FileNotFoundError: If the image file that the PDS3 label names does not
                exist, unless `astrometry` is True.
            OSError: If a file cannot be read.
        """

        GalileoSSI.initialize()  # define everything the first time through; use defaults
                                 # unless initialize() was called explicitly.

        # Default Pds3Label options are method="fast"; the parser's source text is never
        # wanted, in label_dict or anywhere else
        pds3_options = oops.Host._pds3_options_with_defaults(pds3_options)
        pds3_options['sources'] = False

        fileinfo = GalileoSSI._parse_fileinfo(fileinfo, pds3_options=pds3_options)

        # Read the image file unless this is an astrometry read; never read the label
        # beside an image file
        label_dict = None
        vicar = None
        if isinstance(fileinfo, Pds3Label):
            label_dict = _ssi_labels.label_json(fileinfo.as_dict())
            filepath = GalileoSSI._image_path(fileinfo, find=not astrometry)
            if not astrometry:
                vicar = GalileoSSI._read_vicar_image(filepath)
        else:
            filepath = FCPath(fileinfo.filepath)
            vicar = fileinfo
            if not astrometry and not isinstance(vicar, VicarImage):
                vicar = GalileoSSI._read_vicar_image(filepath)

        header_dict = None if vicar is None else _ssi_labels.header_json(vicar)
        data = None if astrometry else vicar.data_2d

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

        snapshot = GalileoSSI._make_snapshot(dict_, filepath=filepath, data=data,
                                             full_fov=full_fov,
                                             return_all_planets=return_all_planets,
                                             overrides=overrides, check_pointing=True,
                                             fetched=vicar is not None)
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
        each later parser of fast, strict, and loose. A VICAR file's header is read
        without its array.

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
            oops.host.PdsHostError: If no parser accepts the PDS3 label.
            oops.host.VicarHostError: If the VICAR header cannot be parsed.
            oops.host.WrongHostError: If the file or the parsed content is not that of a
                Galileo SSI image.
            OSError: If the file cannot be read.
        """

        def read(method):
            return GalileoSSI._read_fileinfo(fileinfo, formats='PV', astrometry=True,
                                             pds3_options=pds3_options
                                             | {'method': method})

        requested = pds3_options['method']
        if requested in _PDS3_METHODS:
            parsed = None
            for method in _PDS3_METHODS[_PDS3_METHODS.index(requested):]:
                try:
                    parsed = read(method)
                    break
                except (oops.host.PdsHostError, oops.OopsValueError) as err:
                    error = err
            if parsed is None:
                raise oops.host.PdsHostError(f'no PDS3 parser can read {fileinfo}: '
                                             f'{error}') from error
        else:
            parsed = read(requested)

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
            FCPath: The file named by the label's IMAGE pointer, in the directory of the
            label. If `find` is True and no file of that name exists, the same name in
            upper or lower case is used if it exists.

        Raises:
            oops.host.PdsHostError: If the label has no IMAGE pointer.
        """

        name = label.get('^IMAGE')
        if not isinstance(name, str):
            raise oops.host.PdsHostError(f'no ^IMAGE pointer in Galileo SSI PDS3 label '
                                         f'{label.filepath}')
        parent = FCPath(label.filepath).parent
        if find:
            for option in dict.fromkeys((name, name.upper(), name.lower())):
                if (parent / option).is_file():
                    return parent / option
        return parent / name

    @staticmethod
    def _read_vicar_image(filepath):
        """The VICAR image in a file, with its header and its full array.

        Parameters:
            filepath (FCPath): The image file.

        Returns:
            VicarImage: The image, whose array is located by its own header.

        Raises:
            oops.host.VicarHostError: If the file is not a valid VICAR image, including
                one whose array is shorter than its header says.
            OSError: If the file cannot be read.
        """

        # A truncated array fails to take the header's shape with a bare ValueError
        try:
            vicar = VicarImage(filepath, strict=False)
        except (VicarError, ValueError) as err:
            raise oops.host.VicarHostError(f'{err}: {filepath}') from err
        return vicar

    @staticmethod
    def from_index(filepath, *, supplement=None, row_dicts=None, pds3_options=None,
                   filter=None, full_fov=False, return_all_planets=False, select=None,
                   target=None, lightsource=None, path=None, frame=None, cadence=None,
                   fov=None, calibrations=None, timeshift=None, navigation=None,
                   tracker=None, parallel=None, **kwargs):
        """A list of Snapshot objects, one for each row of a Galileo SSI index file.

        Parameters:
            filepath (str | pathlib.Path | FCPath): The full path to a Galileo SSI index
                file or its PDS label.
            supplement (str | pathlib.Path | FCPath, optional): The path to the PDS label
                of a supplemental index file, whose columns are merged into each row. Its
                rows must name the same files in the same order as the index; where both
                tables give a column, the index value is kept. If not specified, the
                directory of `filepath` is searched for a file whose name matches except
                that "_index" is replaced by "_supplemental_index" in upper or lower case.
            row_dicts (list[dict[str, Any]], optional): A list of dictionaries, one per
                row of the index. Include if `filepath` and `supplement` have already been
                read. If provided, `supplement` is ignored.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser.
            filter (dict[str, Any], optional): A dictionary of index column names and
                their values or a set or tuple of values. If provided, only rows of the
                index in which each named column equals the value, falls in the set of
                values, or falls in a range defined by a tuple of two values, will be
                included in the returned list.
            full_fov (bool, optional): If True, each field of view covers the full image
                rather than the cutout window.
            return_all_planets (bool, optional): Include kernels for all planets, not just
                the target of the mission phase.
            select (int | slice | str | tuple[int | slice | str, ...]): Accepted for
                compatibility with :meth:`oops.Host.from_index` and ignored, because each
                row describes a single observation.
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
            **kwargs: Additional keyword arguments; they are accepted and ignored.

        Returns:
            list[Snapshot]: One observation per row of the index, including any row with a
            zero exposure, each with subfields `dict`, `host`, `instrument`, `filter`,
            `time_from_sclk`, `label_target`, `observation_id`, `spice_kernels`,
            `filepath`, `basename`, `image_url`, `spice_to_frame`, `spice_frame_name` and
            `spice_frame_id` inserted, `spacecraft_clock_start_count` if the row gives
            it, and `abspath` if the image file is local or already in the local cache.
            `dict` is the row, merged with its row of the supplement, with every value as
            a Python str, int, float, or list; a PDS3 null such as UNK, a numeric -999,
            and the supplement's CUT_OUT_WINDOW of -1, -1, -1, -1 for an image with no
            cutout window are left out. The tables' other null constants are kept,
            because they also mark real values, such as an exposure of 0. No row is
            dropped, and the C kernels are not checked for pointing.
            `filepath`, `basename`, `image_url`, and `abspath` name the image file, as a
            full path or URL: the row's FILE_SPECIFICATION_NAME, with the label's
            extension ".LBL" replaced by ".IMG" in the same case, in the volume directory
            named by its VOLUME_ID. The volume directories are found from the location of
            `filepath`: beside the volume that holds an index in its INDEX directory;
            under "volumes" in place of "metadata" for an index in an RMS Node metadata
            tree (".../metadata/GO_0xxx/GO_0017/GO_0017_index.lbl"); and otherwise in the
            index's own directory.

        Raises:
            FileNotFoundError: If the supplement does not exist.
            oops.host.HostFormatError: If a row's time must come from the spacecraft clock
                count and its FRAME_DURATION is not one of the SSI frame durations.
            oops.host.HostIndexError: If the supplement does not have the same number of
                rows as the index, if a row's VOLUME_ID or FILE_SPECIFICATION_NAME
                differs from the index's, or if `supplement` is not specified and the
                name of `filepath` contains no "_index".
            oops.host.WrongHostError: If `filepath` is not a Galileo SSI index.
            OSError: If `filepath`, or a `supplement` that is specified, does not exist or
                cannot be read.
        """

        GalileoSSI.initialize()

        # Read and filter the index table
        if row_dicts is None:
            supplement = supplement or oops.Host.SUPPLEMENT
            row_dicts = GalileoSSI._read_index_rows(filepath, supplement=supplement,
                                                    validate=_INDEX_KEYS,
                                                    pds3_options=pds3_options)
        row_dicts = oops.Host._filter_index_rows(row_dicts, filter)

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

        # Create the list of Snapshot objects. FILE_SPECIFICATION_NAME names the label;
        # the Snapshot describes the image beside it.
        volumes = GalileoSSI._volumes_directory(filepath)
        snapshots = []
        for row_dict in row_dicts:
            (stem, _, suffix) = str(row_dict['FILE_SPECIFICATION_NAME']).rpartition('.')
            image = (volumes / str(row_dict['VOLUME_ID']).strip()
                     / (stem + '.' + ('img' if suffix.islower() else 'IMG')))
            obs = GalileoSSI._make_snapshot(_ssi_labels.row_items(row_dict),
                                            filepath=image, data=None,
                                            full_fov=full_fov,
                                            return_all_planets=return_all_planets,
                                            overrides=overrides,
                                            default_mode=_INDEX_DEFAULT_MODE,
                                            fetched=False)
            snapshots.append(obs)

        return snapshots

    @staticmethod
    def _volumes_directory(index_path):
        """The directory that holds the volume directories that an index describes.

        Parameters:
            index_path (str | pathlib.Path | FCPath): The index file or its PDS3 label.

        Returns:
            FCPath: The parent of the volume for an index in a volume's INDEX directory;
            for an index in an RMS Node metadata tree, such as
            ".../metadata/GO_0xxx/GO_0017/GO_0017_index.lbl", the matching volume series
            directory, ".../volumes/GO_0xxx"; otherwise the index's own directory.
        """

        directory = FCPath(index_path).parent
        if directory.name.upper() == 'INDEX':
            return directory.parent.parent
        if directory.parent.parent.name == 'metadata':
            return directory.parent.parent.parent / 'volumes' / directory.parent.name
        return directory

    @staticmethod
    def _make_snapshot(dict_, filepath, data=None, *, full_fov=False,
                       return_all_planets=False, overrides, default_mode=None,
                       check_pointing=False, fetched=True):
        """A Snapshot built from the merged metadata of an image or a row of an index.

        Parameters:
            dict_ (dict): The merged metadata of the labels read by :meth:`from_file`, or
                the items of an index row from `_ssi_labels.row_items()`; either holds no
                null value.
            filepath (FCPath | str): The path to the file, or for an index row, the
                VOLUME_ID and FILE_SPECIFICATION_NAME joined by a slash.
            data (numpy.ndarray, optional): The full image array; None to omit the `data`
                subfield.
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
                scan platform's pointing at the midtime of the exposure, unless a frame
                override is given.
            fetched (bool, optional): False if the image file was not read; `abspath` is
                then left out for a remote file that is not already in the local cache.

        Returns:
            Snapshot: The observation, with the subfields listed for :meth:`from_file`.

        Raises:
            oops.host.HostFormatError: If `dict_` gives no exposure duration, no
                telemetry format, or no time (IMAGE_TIME, or the spacecraft clock count
                and frame duration), if the telemetry format is not one of SSI's, or if
                the time must come from the spacecraft clock count and FRAME_DURATION is
                not one of the SSI frame durations.
            oops.OopsSpiceError: If `check_pointing` is True, no frame override is given,
                and the C kernels hold no pointing for the scan platform at the midtime of
                the exposure.
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

        # Apply standard overrides
        oops.Host._apply_overrides(params, filepath, **overrides)

        # Check the pointing at the middle of the exposure the Snapshot will have, which a
        # cadence override can change
        if using_cks and check_pointing:
            GalileoSSI._check_pointing(params['cadence'].midtime, filepath)

        # An astrometry read of a remote label, and an index row, never fetch the image,
        # so its local path would name a file not in the cache unless an earlier read put
        # it there
        if (not fetched and not isinstance(filepath, str) and not filepath.is_local()
                and not params['abspath'].exists()):
            del params['abspath']

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
    def _check_pointing(time, filepath):
        """Raise OopsSpiceError unless the C kernels give the scan platform's pointing.

        Parameters:
            time (float): The time in seconds TDB.
            filepath (FCPath): The image file, named in the error.

        Raises:
            oops.OopsSpiceError: If the scan platform frame cannot be evaluated at `time`.
                Its message is "no Galileo SSI pointing at <time> for <file>: <code>",
                with the SPICE error code, such as "SPICE(CKINSUFFDATA)", and its cause
                is the exception that cspyce raised.
        """

        frame = oops.Frame.as_frame(_SPICE_FRAME_NAME)
        try:
            frame.transform_at_time(time)
        except _SPICE_ERRORS as err:
            raise oops.OopsSpiceError.from_spice_error(err, description=_NO_POINTING,
                                                       time=time,
                                                       filepath=filepath) from err

    @staticmethod
    def _tstart_from_sclk(sclk, frame_duration, exposure_type, filepath):
        """The start of the exposure derived from the spacecraft clock count.

        The count names the frame in which the image was read out. The shutter opens a
        fixed delay after the start of that frame, which depends on FRAME_DURATION, or
        one frame earlier for an EXTENDED exposure taken before the Jupiter tour's
        phase-2 flight software. The result agrees with the start of the exposure implied
        by IMAGE_TIME to within about 1/3 s for most images.

        Parameters:
            sclk (str): The SPACECRAFT_CLOCK_START_COUNT, as "RRRRRRRR.MM".
            frame_duration (float | str): The FRAME_DURATION in seconds.
            exposure_type (str | None): The EXPOSURE_TYPE; None is taken as "NORMAL".
            filepath (FCPath | str): The image file, named in an error.

        Returns:
            float: The start of the exposure in seconds TDB.

        Raises:
            oops.host.HostFormatError: If FRAME_DURATION is not one of the SSI frame
                durations.
        """

        frame_duration = float(frame_duration)
        nominal = min(_SHUTTER_DELAY, key=lambda key: abs(key - frame_duration))
        if abs(nominal - frame_duration) > _FRAME_DURATION_TOLERANCE:
            raise oops.host.HostFormatError(f'unrecognized Galileo SSI FRAME_DURATION '
                                            f'{frame_duration} for {filepath}')

        tstart = _Galileo.tdb_from_sclk(sclk) + _SHUTTER_DELAY[nominal]
        phase1 = int(sclk.partition('.')[0]) < _PHASE2_FIRST_RIM
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

        A Galileo SSI index label names the host and instrument at its top level, so
        unlike the other detectors, this one never needs to look at the rows.
        """

        return (label.get('INSTRUMENT_HOST_NAME', '').startswith('GALILEO')
                and label.get('INSTRUMENT_ID', '') == 'SSI')

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

        Fills in key information about the camera. Must be called first. After the first
        call, later calls to this function are ignored.

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
