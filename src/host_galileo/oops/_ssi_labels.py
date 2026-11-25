##########################################################################################
# host_galileo/oops/_ssi_labels.py
##########################################################################################
"""The metadata of an SSI image, from its PDS3 label, its VICAR header, or both.

Each source is first reduced to its metadata items. The PDS3 label gives its values under
their own names. The VICAR header gives each value that the host converts into a PDS3
item only under the PDS3 name, with the type the PDS3 label gives it, and every other
value under its own name. A value that a source gives as "no value" is not an item of
that source. Where both sources give an item, :func:`merge_items` keeps the value of the
item's one fixed source.
"""

import datetime as dt

from oops.host import is_pds3_null

# The value of a numeric VICAR keyword that the header leaves unset, such as SMRAZ
_VICAR_FILL = -999

# The year of the VICAR time fields of a frame that has no time, as in 48 of the 77
# RAW_CAL frames of GO_0002 and GO_0003 whose labels give IMAGE_TIME = UNK
_VICAR_NO_YEAR = -32768

# The orbits of the Jupiter tour, from J0 to orbit 35. The labels of GO_0003-GO_0009 and
# GO_0016 give the ORBIT_NUMBER of every frame as two characters of its observation ID,
# byte-swapped, read as a 16-bit integer, such as 12613 for "E1" ("1E" = 0x3145); 5,213
# of the 14,672 products. These frames precede the Jupiter tour and have no orbit number,
# and the DATASET.CAT of GO_0016 lists ORBIT_NUMBER as not applicable. A value outside the
# tour is therefore treated as missing; this rule is listed with the other label
# differences in SETI/rms-host-galileo#41.
_ORBITS = range(36)

# The label items that the PDS3 label writes as a set although their order carries
# meaning; each is taken from the parser's list, which keeps the label's order
_ORDERED_SETS = ('CUT_OUT_WINDOW', 'TRUTH_WINDOW')

# The PDS3 keywords and the VICAR keywords that give the same whole value. Some values
# have the same name in both, such as CUT_OUT_WINDOW and TARGET_CENTER_DISTANCE; they need
# no entry. The phase-1 headers of GO_0002-GO_0016 spell PARTITION as PARTITIO.
_VICAR_PDS3_PAIRS = [
    ('SPACECRAFT_CLOCK_CNT_PARTITION', 'PARTITION'                  ),
    ('SPACECRAFT_CLOCK_CNT_PARTITION', 'PARTITIO'                   ),
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

# The FILTER_NAME of the PDS3 label, indexed by its FILTER_NUMBER, the FILTER of the VICAR
# header
_FILTER_NAMES = ['CLEAR', 'GREEN', 'RED', 'VIOLET', 'IR-7560', 'IR-9680', 'IR-7270',
                 'IR-8890']

# The FRAME_DURATION of the PDS3 label in seconds, keyed by the RATE of the VICAR header
_FRAME_DURATIONS = {1: 2.333, 2: 8.667, 3: 30.333, 4: 60.667, 5: 15.166}

# The GAIN_MODE_ID of the PDS3 label, keyed by the GAIN of the VICAR header
_GAIN_MODES = {1: '400K', 2: '100K', 3: '40K', 4: '10K'}

# The VICAR time fields, after their prefix, and the PDS3 item each set of them gives
_TIME_FIELDS = ('YEAR', 'DAY', 'HOUR', 'MIN', 'SEC', 'MSEC')
_TIME_ITEMS = (('IMAGE_TIME', 'SCET'), ('EARTH_RECEIVED_TIME', 'ERT'))

# The resolution rules: for each item that the two sources give differently, the one
# source whose value is kept for every image. An item without a rule keeps the PDS3
# label's value. A value missing from the winning source is not a difference; the other
# source's value is used. The counts come from comparing the two labels of all 14,672 SSI
# products of GO_0002-GO_0023 ("of" counts the products that give the item in both), and
# the SPICE checks evaluate the host's own kernels at each image's IMAGE_TIME. Every rule,
# with its evidence and the images on which its source is wrong, is recorded in
# SETI/rms-host-galileo#41.
#
# MIPL wrote the VICAR header first, updating its pointing and geometry from the latest
# SPICE kernels (CATLABEL), and generated the PDS3 labels from it afterwards (CDGEN), as
# CATALOG/DATASET.CAT describes. Where the two differ only because the label rounded the
# header's number to fewer digits, the header holds the value as computed.
_LABEL = 'label'
_HEADER = 'header'
_RULES = {
    # The label is the header rounded to fewer digits in every product that differs;
    # never beyond rounding.
    'CONE_ANGLE'                    : _HEADER,  # 2,543 of 14,422
    'DECLINATION'                   : _HEADER,  # 13,642 of 14,573
    'EMISSION_ANGLE'                : _HEADER,  # 9,573 of 10,545
    'ENTROPY'                       : _HEADER,  # 14,406 of 14,672
    'INCIDENCE_ANGLE'               : _HEADER,  # 8,864 of 10,545
    'INTERCEPT_POINT_LATITUDE'      : _HEADER,  # 3,660 of 3,853
    'INTERCEPT_POINT_LINE'          : _HEADER,  # 1,084 of 3,753
    'INTERCEPT_POINT_LINE_SAMPLE'   : _HEADER,  # 1,066 of 3,753
    'INTERCEPT_POINT_LONGITUDE'     : _HEADER,  # 848 of 3,902
    'LOCAL_HOUR_ANGLE'              : _HEADER,  # 419 of 3,847
    'NORTH_AZIMUTH'                 : _HEADER,  # 6,948 of 11,360
    'PHASE_ANGLE'                   : _HEADER,  # 8,720 of 10,545
    'SUB_SOLAR_AZIMUTH'             : _HEADER,  # 213 of 11,360
    'SUB_SPACECRAFT_AZIMUTH'        : _HEADER,  # 4,898 of 11,360
    'SUB_SPACECRAFT_LONGITUDE'      : _HEADER,  # 778 of 3,753
    'TWIST_ANGLE'                   : _HEADER,  # 4,260 of 14,573

    # Rounding in 5,115 of 14,573; in the other one, GO_0007/REDO/C0059466445R, the
    # header writes the same direction as -45.4753, outside the documented 0-360.
    'RIGHT_ASCENSION'               : _HEADER,

    # The SSI exposures are fixed commanded steps such as 4 1/6 and 45 5/6 ms. The
    # header gives them to within 1 us (4.167, 45.8333; 12.5003 for 12.5 in 388), the
    # label to two decimal places, up to 5 us off (4.17); 7,576 of 14,672 differ.
    'EXPOSURE_DURATION'             : _HEADER,

    # Checked against SPICE: where the two differ beyond rounding, the header matches and
    # the label is wrong, mostly in GO_0017-GO_0020 SML_SATS frames, by up to 200,000 km
    # (81 products) or 108 degrees (73 for SUB_SOLAR_LONGITUDE, 39 for
    # SUB_SOLAR_LATITUDE). The label's SUB_SPACECRAFT_LINE and _SAMPLE are double the
    # header's in summed (HIS) frames; the header matches the projected target center
    # in 86 of the 109 frames that differ beyond rounding and settle it. The 11
    # SUB_SPACECRAFT_LATITUDE values of GO_0020 E11 that differ by 0.0015 degrees match
    # SPICE equally well.
    'SUB_SOLAR_LATITUDE'            : _HEADER,  # 41 of 3,880 beyond rounding
    'SUB_SOLAR_LONGITUDE'           : _HEADER,  # 90 of 3,880
    'SUB_SPACECRAFT_LATITUDE'       : _HEADER,  # 11 of 3,753
    'SUB_SPACECRAFT_LINE'           : _HEADER,  # 134 of 3,753
    'SUB_SPACECRAFT_LINE_SAMPLE'    : _HEADER,  # 148 of 3,753
    'TARGET_CENTER_DISTANCE'        : _HEADER,  # 139 of 3,880

    # DATASET.CAT: the label's compression ratios of an on-chip mosaic's child images are
    # copied from the parent image, as rough approximations; the header's are the
    # image's own.
    'ENCODING_COMPRESSION_RATIO'    : _HEADER,  # 303 of 3,592 beyond rounding
    'ENCODING_MAX_COMPRESSION_RATIO': _HEADER,  # 178 of 3,592
    'ENCODING_MIN_COMPRESSION_RATIO': _HEADER,  # 247 of 3,592

    # Checked against SPICE: the label is nearer in 5,419 of 6,560 differing Earth-flyby
    # products and right by definition in 301 Moon images (the planet the Moon orbits is
    # the Earth; the header gives the Moon's own distance); the header is nearer in 1,737
    # of 1,748 tour products, where the label is off by a median 64 km. 8,745 of 11,900
    # differ beyond rounding.
    'CENTRAL_BODY_DISTANCE'         : _LABEL,

    # The header says BARC INFORMATION PRESERVING where the label says BARC RATE CONTROL
    # (74 of 4,003); 56 of the 74 have truncated samples, which only rate control
    # produces.
    'COMPRESSION_TYPE'              : _LABEL,

    # The labels correct the header: 150 respell BLACK_SKY or NON_SCIENCE with a space,
    # and 51 name a ring image's target J RINGS where the header names JUPITER (each with
    # the label's RING_RADIUS and the header's RAD of 180,000 km); 201 of 14,672.
    'TARGET_NAME'                   : _LABEL,
}


def label_items(label_dict):
    """The metadata items of a PDS3 label.

    Parameters:
        label_dict (dict): The label, as given by ``Pds3Label.as_dict()``.

    Returns:
        dict: Every key of `label_dict` with its value, except a key whose value is None
        or a PDS3 null text (as ``oops.host.is_pds3_null()`` decides), and an
        ORBIT_NUMBER outside 0-35, the orbits of the Jupiter tour. CUT_OUT_WINDOW and
        TRUTH_WINDOW, which the label writes as sets, are given as lists in the label's
        order.
    """

    items = {key: value for (key, value) in label_dict.items()
             if value is not None and not is_pds3_null(value)}
    for key in _ORDERED_SETS:
        if key in items:
            items[key] = items[key + '_list']
    if items.get('ORBIT_NUMBER') not in _ORBITS:
        items.pop('ORBIT_NUMBER', None)
    return items


def header_items(header):
    """The metadata items of a VICAR header, under PDS3 names where converted.

    Parameters:
        header (VicarLabel): The VICAR header.

    Returns:
        dict: Every parameter of the header, keyed by its name, except a value of -999 or
        a PDS3 null text; a repeated parameter has the value of its first occurrence. A
        value that the host converts into a PDS3 item is given only under the PDS3 name:
        the keywords that name the same value (such as EXP for EXPOSURE_DURATION, FILTER
        for FILTER_NUMBER, and PA for OBSERVATION_ID); RIM and MOD91 as
        SPACECRAFT_CLOCK_START_COUNT, written "RRRRRRRR.MM"; the SCET and ERT fields as
        IMAGE_TIME and EARTH_RECEIVED_TIME, each a ``datetime.datetime``, or neither if
        the year is -32768; and the codes RATE and GAIN as FRAME_DURATION and
        GAIN_MODE_ID. FILTER_NAME comes from the filter number, and EXPOSURE_TYPE from
        the last flag of MOFIBE, or of FIBE if there is no MOFIBE. MOD10, MOD8, FIBE, and
        MOFIBE keep their own names, since no PDS3 item holds them.
    """

    items = {key: value for (key, value) in header.as_dict(unique=False).items()
             if value != _VICAR_FILL and not is_pds3_null(value)}

    for (pds3_key, vicar_key) in _VICAR_PDS3_PAIRS:
        if vicar_key in items:
            items[pds3_key] = items.pop(vicar_key)

    if 'RIM' in items and 'MOD91' in items:
        items['SPACECRAFT_CLOCK_START_COUNT'] = (f'{items.pop("RIM"):08d}.'
                                                 f'{items.pop("MOD91"):02d}')

    _set_coded(items, 'FILTER_NAME', 'FILTER_NUMBER', dict(enumerate(_FILTER_NAMES)))
    _set_coded(items, 'FRAME_DURATION', 'RATE', _FRAME_DURATIONS, pop=True)
    _set_coded(items, 'GAIN_MODE_ID', 'GAIN', _GAIN_MODES, pop=True)

    # The last flag of FIBE, or of MOFIBE in the phase-2 headers of the Jupiter tour,
    # marks an extended exposure
    flags = items.get('MOFIBE', items.get('FIBE'))
    if flags:
        items['EXPOSURE_TYPE'] = 'EXTENDED' if flags[-1] == '1' else 'NORMAL'

    for (pds3_key, prefix) in _TIME_ITEMS:
        time = _vicar_time([items.pop(prefix + name, None) for name in _TIME_FIELDS])
        if time is not None:
            items[pds3_key] = time

    return items


def _set_coded(items, pds3_key, vicar_key, names, *, pop=False):
    """Give the PDS3 item that a VICAR code stands for, if the code is known.

    Parameters:
        items (dict): The header's items, modified in place.
        pds3_key (str): The name of the PDS3 item.
        vicar_key (str): The name of the item that holds the code.
        names (dict): The PDS3 value for each known code.
        pop (bool, optional): True to remove the code once it is converted.
    """

    code = items.get(vicar_key)
    if code in names:
        items[pds3_key] = names[code]
        if pop:
            del items[vicar_key]


def _vicar_time(fields):
    """The time given by a set of VICAR time fields.

    Parameters:
        fields (list[int | None]): The year, day of the year, hour, minute, second, and
            millisecond, each None if the header leaves it out.

    Returns:
        datetime.datetime | None: The time, or None if the header leaves out a field or
        gives the year -32768.
    """

    if None in fields or fields[0] == _VICAR_NO_YEAR:
        return None
    (year, doy, hour, minute, second, msec) = fields
    return (dt.datetime(year, 1, 1, hour, minute, second, msec * 1000)
            + dt.timedelta(days=doy - 1))


def merge_items(label, header):
    """The merged metadata of the sources that were read.

    Parameters:
        label (dict | None): The PDS3 label's items from :func:`label_items`, or None if
            the label was not read.
        header (dict | None): The VICAR header's items from :func:`header_items`, or None
            if the header was not read.

    Returns:
        dict: Every item of either source. An item that both give with different values
        has the value of the source that its rule names, or the PDS3 label's value if it
        has no rule.
    """

    merged = dict(header or {})
    merged.update(label or {})
    if label and header:
        for key in label.keys() & header.keys():
            if _RULES.get(key) == _HEADER:
                merged[key] = header[key]
    return merged

##########################################################################################
