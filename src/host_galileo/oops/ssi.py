##########################################################################################
# host_galileo/oops/ssi.py
##########################################################################################

import datetime as dt
import re

import julian
import cspyce

from filecache import FCPath
from pdsparser import Pds3Label
from vicar     import VicarImage

import oops
from ._galileo import _Galileo

__all__ = ['GalileoSSI']

# The SSI camera frame is the scan platform frame, whose C kernels are sampled at discrete
# spacecraft clock ticks. SSI images are spaced as closely as 1 unit in the file name,
# which corresponds to 80 clock ticks. Therefore, we use a tolerance of +/-40.
_SPICE_FRAME_NAME = 'GLL_SCAN_PLATFORM'
_SPICE_FRAME_ID = -77001
_TICK_TOLERANCE = 40

# The NAIF ID of the SSI instrument, as used in its instrument kernel
_NAIF_ID = -77036

# The CUT_OUT_WINDOW value written in the supplemental index when the label has none
_NO_WINDOW = [-1, -1, -1, -1]

# The columns that identify the same image in an index and its supplement. Both tables
# also give SPACECRAFT_CLOCK_START_COUNT, but the supplement writes it in the SPICE form
# "RRRRRRRR:MM:S:T" and the index in the label form "RRRRRRRR.MM", so the default
# validation of oops.Host._read_index_rows(), which compares every shared column, rejects
# every real pair.
_INDEX_KEYS = ['VOLUME_ID', 'FILE_SPECIFICATION_NAME']

# The dimensions of an SSI image in pixels, as the instrument kernel must give them
_SSI_PIXELS = 800

# Seconds from the start of the frame named by SPACECRAFT_CLOCK_START_COUNT to the opening
# of the shutter, keyed by FRAME_DURATION in seconds. These are the median offsets between
# the start of the exposure (IMAGE_TIME, the middle of the exposure, minus half of
# EXPOSURE_DURATION) and the clock count across the 14,585 images of GO_0002-GO_0023 that
# give both; about 90% of images agree to within 1/3 s, the resolution of a count given to
# the mod-91 field. An EXTENDED exposure opens its shutter one frame earlier, in the frame
# before the one the count names.
_SHUTTER_DELAY = {
     2.333: 0.225,
     8.667: 1.149,
    15.167: 1.150,
    30.333: 2.817,
    60.667: 6.483,
}
_FRAME_DURATION_TOLERANCE = 0.05

# Target name as found in PDS label -> correct name
_TARGET_NAME_REPAIRS = {
    'BLACK_SKY'  : 'NONE',
    'NON_SCIENCE': 'NONE',
    'STAR'       : 'NONE',
    'HD151288'   : 'HD 151288',
    'J RINGS'    : 'JUPITER_RING_PLANE',
}

# These are star targets and need to be defined as LightSources
# "BET_ARI", "DEL_PSC", "GAM_ORI", "HD 151288", "PLEIADES", "SIRIUS", "TAU_CET"
_TARGET_STARS = re.compile(r'(HD 151288|PLEIADES|SIRIUS|[A-Z]{3}_[A-Z]{3})$')


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

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | Pds3Label | VicarLabel | VicarImage):
                The path to the PDS3 label or the VICAR image file, or its parsed content
                as a Pds3Label, VicarLabel, or VicarImage. The metadata of a VICAR file
                are read from its own header; no PDS3 label is needed.
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
            Snapshot: The observation, with subfields `dict`, `host`, `instrument`,
            `filter`, `time_from_sclk`, `label_target`, `observation_id`,
            `spice_kernels`, `filepath`, `basename`, `spice_to_frame`, `spice_frame_name`,
            `spice_frame_id`, `abspath` and `image_url` inserted, and `data` unless
            `astrometry` is True. `dict` is the PDS3 label, or the VICAR header with its
            values also given under their PDS3 keyword names, and `time_from_sclk` is True
            if the label gives no IMAGE_TIME, so that the time comes from the spacecraft
            clock count.

        Raises:
            oops.host.PdsHostError: If the PDS3 label cannot be parsed.
            oops.host.VicarHostError: If the VICAR header cannot be parsed.
            oops.host.HostFormatError: If the time must come from the spacecraft clock
                count and FRAME_DURATION is not one of the SSI frame durations.
            OSError: If the file cannot be read.
        """

        GalileoSSI.initialize()  # define everything the first time through; use defaults
                                 # unless initialize() was called explicitly.

        # Default Pds3Label options are method="fast", sources=True
        pds3_options = oops.Host._pds3_options_with_defaults(pds3_options)
        pds3_options['sources'] = True      # required

        fileinfo = GalileoSSI._read_fileinfo(fileinfo, formats='PV',
                                             astrometry=astrometry,
                                             pds3_options=pds3_options)

        # Read the data if necessary; harmonize the label content
        data = None
        if isinstance(fileinfo, Pds3Label):
            if not astrometry:
                data = oops.host.read_pds3_image_array(fileinfo)
            dict_ = fileinfo.as_dict()
            # Fix the lists incorrectly labeled as sets
            for key in ('CUT_OUT_WINDOW', 'TRUTH_WINDOW'):
                if key in dict_:
                    dict_[key] = dict_[key + '_list']
        else:
            if not astrometry:
                if not isinstance(fileinfo, VicarImage):
                    # A VicarImage built from a VicarLabel has no array; re-read the file
                    fileinfo = VicarImage(fileinfo.filepath, strict=False)
                data = fileinfo.data_2d
            dict_ = GalileoSSI._vicar_to_pds3(fileinfo)

        filepath = FCPath(fileinfo.filepath)

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
                                             overrides=overrides)
        return snapshot

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
            `filepath`, `basename`, `spice_to_frame`, `spice_frame_name` and
            `spice_frame_id` inserted. `dict` is the row, merged with its row of the
            supplement, and `filepath` is the VOLUME_ID and FILE_SPECIFICATION_NAME
            joined by a slash.

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

        # Create the list of Snapshot objects
        snapshots = []
        for row_dict in row_dicts:
            fpath = row_dict['VOLUME_ID'] + '/' + row_dict['FILE_SPECIFICATION_NAME']
            obs = GalileoSSI._make_snapshot(row_dict, filepath=fpath, data=None,
                                            full_fov=full_fov,
                                            return_all_planets=return_all_planets,
                                            overrides=overrides)
            snapshots.append(obs)

        return snapshots

    @staticmethod
    def _make_snapshot(dict_, filepath, data=None, *, full_fov=False,
                       return_all_planets=False, overrides):
        """A Snapshot built from a PDS3 label or a row of an index.

        Parameters:
            dict_ (dict): The PDS3 label as a dictionary, or an index row.
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

        Returns:
            Snapshot: The observation, with the subfields listed for :meth:`from_file`.

        Raises:
            oops.host.HostFormatError: If the time must come from the spacecraft clock
                count and FRAME_DURATION is not one of the SSI frame durations.
        """

        # Note that PDS3 labels and index rows use the same names.
        texp = dict_['EXPOSURE_DURATION'] / 1000.

        # IMAGE_TIME is the middle of the exposure. Some RAW_CAL frames in GO_0002 and
        # GO_0003 have no IMAGE_TIME, but every label carries the spacecraft clock count.
        # Never fall back to a placeholder time: it silently places the frame at the
        # J2000 epoch.
        if dict_['IMAGE_TIME'] == 'UNK':
            tstart = GalileoSSI._tstart_from_sclk(dict_)
            time_from_sclk = True
        else:
            tmid = julian.tdb_from_tai(julian.tai_from_iso(dict_['IMAGE_TIME']))
            tstart = tmid - texp / 2.
            time_from_sclk = False

        mode = dict_.get('TELEMETRY_FORMAT_ID', 'NONE')
        filter_ = dict_['FILTER_NAME']

        label_target = dict_.get('TARGET_NAME')
        label_target = _TARGET_NAME_REPAIRS.get(label_target, label_target)
        if _TARGET_STARS.match(label_target):
            label_lightsource = oops.lightsource.star_lookup(label_target)
            label_target = 'NONE'
        else:
            label_lightsource = 'SUN'

        # Make sure the SPICE frame is defined unless a custom frame was given
        using_cks = overrides['frame'] is None
        if using_cks:
            GalileoSSI._define_camera_frame()

        # The cutout window is (line, sample, lines, samples), with a one-based origin
        window = dict_.get('CUT_OUT_WINDOW')
        if window is not None and list(window) == _NO_WINDOW:
            window = None

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
            'target'      : label_target,
            'lightsource' : label_lightsource,
            'calibrations': [],         # TBD
        }

        # Apply standard overrides
        oops.Host._apply_overrides(params, filepath, **overrides)

        obs = oops.observation.Snapshot(('v', 'u'), **params)
        if data is not None:
            obs.insert_subfield('data', data)

        # Extra attributes: use standardized names wherever practical
        obs.insert_subfield('dict', dict_)
        obs.insert_subfield('host', 'Galileo')
        obs.insert_subfield('instrument', 'SSI')
        obs.insert_subfield('filter', filter_)
        obs.insert_subfield('time_from_sclk', time_from_sclk)
        obs.insert_subfield('label_target', dict_['TARGET_NAME'])
        obs.insert_subfield('observation_id', dict_['OBSERVATION_ID'])
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
    def _tstart_from_sclk(dict_):
        """The start of the exposure derived from the spacecraft clock count.

        The count names the frame in which the image was read out. The shutter opens a
        fixed delay after the start of that frame, which depends on FRAME_DURATION, or
        one frame earlier for an EXTENDED exposure. The result agrees with the start of
        the exposure implied by IMAGE_TIME to within about 1/3 s for most images.

        Parameters:
            dict_ (dict): The PDS3 label or index row. It must give
                SPACECRAFT_CLOCK_START_COUNT and FRAME_DURATION; EXPOSURE_TYPE, if absent,
                is taken as "NORMAL".

        Returns:
            float: The start of the exposure in seconds TDB.

        Raises:
            oops.host.HostFormatError: If FRAME_DURATION is not one of the SSI frame
                durations.
        """

        frame_duration = float(dict_['FRAME_DURATION'])
        nominal = min(_SHUTTER_DELAY, key=lambda key: abs(key - frame_duration))
        if abs(nominal - frame_duration) > _FRAME_DURATION_TOLERANCE:
            raise oops.host.HostFormatError(f'unrecognized Galileo SSI FRAME_DURATION: '
                                            f'{frame_duration}')

        tstart = (_Galileo.tdb_from_sclk(dict_['SPACECRAFT_CLOCK_START_COUNT'])
                  + _SHUTTER_DELAY[nominal])
        if dict_.get('EXPOSURE_TYPE', 'NORMAL') == 'EXTENDED':
            tstart -= nominal

        return tstart

    # A table of VICAR keywords and their PDS3 label equivalents
    # These VICAR keywords are the same in the PDS3 labels:
    # TARGET_CENTER_DISTANCE
    # ICT_DESPIKE_THRESHOLD
    # CUT_OUT_WINDOW
    # HUFFMAN_TABLE_TYPE
    # SUB_SOLAR_LATITUDE
    # SUB_SOLAR_LONGITUDE
    # SUB_SPACECRAFT_LATITUDE
    # SUB_SPACECRAFT_LONGITUDE
    # SUB_SPACECRAFT_LINE
    # INTERCEPT_POINT_LINE
    # INTERCEPT_POINT_LINE_SAMPLE
    # ENTROPY
    _VICAR_PDS3_PAIRS = [
        ('SPACECRAFT_CLOCK_CNT_PARTITION', 'PARTITION'                  ),
        ('IMAGE_ID'                      , 'PICNO'                      ),
        ('OBSERVATION_ID'                , 'PA'                         ),
        ('TARGET_NAME'                   , 'TARGET'                     ),
        ('NTV_TIME_FROM_CLOSEST_APPROACH', 'TCA'                        ),
        ('FILTER_NUMBER'                 , 'FILTER'                     ),
        ('EXPOSURE_DURATION'             , 'EXP'                        ),
        ('COMPRESSION_TYPE'              , 'ENCODING_TYPE'              ),
        ('TELEMETRY_FORMAT_ID'           , 'TLMFMT'                     ),
        ('CENTRAL_BODY_DISTANCE'         , 'PLRANGE'                    ),
        ('SUB_SPACECRAFT_LINE_SAMPLE'    , 'SUB_SPACECRAFT_SAMPLE'      ),
        ('SUB_SPACECRAFT_AZIMUTH'        , 'SCAZ'                       ),
        ('TWIST_ANGLE'                   , 'TWIST'                      ),
        ('CONE_ANGLE'                    , 'CONE'                       ),
        ('RIGHT_ASCENSION'               , 'RA'                         ),
        ('DECLINATION'                   , 'DEC'                        ),
        ('NORTH_AZIMUTH'                 , 'NORAZ'                      ),
        ('SMEAR_AZIMUTH'                 , 'SMRAZ'                      ),
        ('SUB_SOLAR_AZIMUTH'             , 'SUNAZ'                      ),
        ('INCIDENCE_ANGLE'               , 'INA'                        ),
        ('EMISSION_ANGLE'                , 'EMA'                        ),
        ('LOCAL_HOUR_ANGLE'              , 'HRA'                        ),
        ('INTERCEPT_POINT_LATITUDE'      , 'LAT'                        ),
        ('INTERCEPT_POINT_LONGITUDE'     , 'LON'                        ),
        ('ENCODING_COMPRESSION_RATIO'    , 'COMPRESSION_RATIO'          ),
        ('ENCODING_MIN_COMPRESSION_RATIO', 'MINIMUM_COMPRESSION_RATIO'  ),
        ('ENCODING_MAX_COMPRESSION_RATIO', 'MAXIMUM_COMPRESSION_RATIO'  ),
        ('ICT_ZIGZAG_PATTERN'            , 'ZIGZAG_PATTERN'             ),
        ('ICT_QUANTIZATION_STEP_SIZE'    , 'QUANTIZATION_STEP_SIZE'     ),
        ('CMPRS_QUANTZ_TBL_ID'           , 'QUANTIZATION_MATRIX_NAME'   ),
        ('PHASE_ANGLE'                   , 'PHA'                        ),
    ]

    _FILTER_NAMES = ['CLEAR', 'GREEN', 'RED', 'VIOLET', 'IR-7560', 'IR-9680', 'IR-7270',
                     'IR-8890']

    @staticmethod
    def _vicar_to_pds3(vic):
        """A VICAR label converted to match the content of a PDS3 label."""

        dict_ = vic.as_dict(unique=False)
        dict_['SPACECRAFT_NAME'] = 'GALILEO ORBITER'
        dict_['INSTRUMENT_NAME'] = 'SOLID STATE IMAGING SYSTEM'
        for (pds3_key, vicar_key) in GalileoSSI._VICAR_PDS3_PAIRS:
            vicar_value = dict_.get(vicar_key, -999.)
            dict_[pds3_key] = 'UNK' if vicar_value == -999. else vicar_value

        sclk = str(dict_['RIM']).rjust(8, '0') + '.' + str(dict_['MOD91']).rjust(2, '0')
        dict_['SPACECRAFT_CLOCK_START_COUNT'] = sclk

        dict_['IMAGE_TIME'] = GalileoSSI._vicar_dict_to_datetime(dict_, 'SCET')
        if 'ERTYEAR' in dict_:
            dict_['EARTH_RECEIVED_TIME'] = GalileoSSI._vicar_dict_to_datetime(dict_,
                                                                              'ERT')

        orbit = dict_['OBSERVATION_ID'][:2]
        dict_['ORBIT_NUMBER'] = int(orbit if orbit.isdigit() else orbit[1])

        dict_['FILTER_NAME'] = GalileoSSI._FILTER_NAMES[dict_['FILTER_NUMBER']]
        return dict_

    @staticmethod
    def _vicar_dict_to_datetime(dict_, prefix):
        """Convert a set of time fields in a VICAR header to a datetime."""

        year = dict_[prefix + 'YEAR']
        doy = dict_[prefix + 'DAY']
        h = dict_[prefix + 'HOUR']
        m = dict_[prefix + 'MIN']
        s = dict_[prefix + 'SEC']
        ms = dict_[prefix + 'MSEC']

        return dt.datetime(year, 1, 1, h, m, s, ms*1000) + dt.timedelta(days=(doy-1))

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
        GalileoSSI._INSTRUMENT_KERNEL = _Galileo.spice_instrument_kernel('SSI')[0]

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
