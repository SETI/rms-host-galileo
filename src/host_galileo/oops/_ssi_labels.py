##########################################################################################
# host_galileo/oops/_ssi_labels.py
##########################################################################################
"""The metadata of an SSI image, from its PDS3 label, its VICAR header, or both.

Each source is first reduced to its metadata items: the PDS3 label's values, and the
VICAR header's values together with the same values under their PDS3 names, each spelled
as the label spells it. A value that a source gives as "no value" (a PDS3 null such as
UNK, or a VICAR fill such as PICNO "?" or -999) is not an item of that source. Where both
sources give an item, :func:`merge_items` keeps one value by the rule for that item.
"""

import datetime as dt
import re

# The PDS3 constants that stand for a missing value
_PDS3_NULLS = ('UNK', 'N/A', 'NULL')

# The value of a numeric VICAR keyword that the header leaves unset, such as SMRAZ
_VICAR_FILL = -999.

# The PICNO of a VICAR header that gives no picture number
_NO_PICNO = '?'

# The PDS3 keywords and the VICAR keywords that give the same value. Some values have
# the same name in both, such as CUT_OUT_WINDOW and TARGET_CENTER_DISTANCE; they need no
# entry. The phase-1 headers of GO_0002-GO_0016 spell PARTITION as PARTITIO.
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

# The FILTER_NAME of the PDS3 label, indexed by the FILTER of the VICAR header
_FILTER_NAMES = ['CLEAR', 'GREEN', 'RED', 'VIOLET', 'IR-7560', 'IR-9680', 'IR-7270',
                 'IR-8890']

# The FRAME_DURATION of the PDS3 label in seconds, keyed by the RATE of the VICAR header
_FRAME_DURATIONS = {1: 2.333, 2: 8.667, 3: 30.333, 4: 60.667, 5: 15.166}

# The GAIN_MODE_ID of the PDS3 label, keyed by the GAIN of the VICAR header
_GAIN_MODES = {1: '400K', 2: '100K', 3: '40K', 4: '10K'}

# The orbits of the Jupiter tour, as the code that begins an observation ID: the orbit
# number as two digits from orbit 10 on, and before that the initial of the orbit's
# encounter target followed by its single digit. The IDs of the earlier mission phases,
# including those of the Earth flybys, which begin "E1" and "E2", name no orbit.
_ORBIT_CODE = re.compile(r'\d\d|J0|G1|G2|C3|E4|J5|E6|G7|G8|C9')

# The orbit numbers of the Jupiter tour, from J0 to orbit 35
_ORBITS = range(36)

# The resolution rules, for the items that the two sources give differently. Each item
# names the source whose value is kept, and an item without a rule keeps the PDS3 label's
# value, the archive's curated description of the product. Before a rule applies, a
# number that is the other source's number rounded gives way to it (see _is_rounded). The
# counts come from comparing the two labels of all 14,672 SSI products of GO_0002-GO_0023;
# each is the number of products that give the item in both and differ beyond rounding,
# out of the number that give it in both. The data give no reason to prefer the header
# for any of these, so every rule keeps the label. Every difference, and every item that
# only one source gives, is recorded in #41.
#
# These items differ only by rounding, so the more precise value is always kept, with the
# number of products that differ: CONE_ANGLE 2,543, DECLINATION 13,642, EMISSION_ANGLE
# 9,573, ENTROPY 14,406, EXPOSURE_DURATION 7,576 (the label gives two decimal places, the
# header up to four), INCIDENCE_ANGLE 8,864, INTERCEPT_POINT_LATITUDE 3,660,
# INTERCEPT_POINT_LINE 1,084, INTERCEPT_POINT_LINE_SAMPLE 1,066, INTERCEPT_POINT_LONGITUDE
# 848, LOCAL_HOUR_ANGLE 419, NORTH_AZIMUTH 6,948, PHASE_ANGLE 8,720, SUB_SOLAR_AZIMUTH
# 213, SUB_SPACECRAFT_AZIMUTH 4,898, SUB_SPACECRAFT_LONGITUDE 778, and TWIST_ANGLE 4,260.
_LABEL = 'label'
_HEADER = 'header'
_RULES = {
    'CENTRAL_BODY_DISTANCE'         : _LABEL,   # 8,745 of 12,648
    'COMPRESSION_TYPE'              : _LABEL,   # 74 of 4,003
    'ENCODING_COMPRESSION_RATIO'    : _LABEL,   # 303 of 3,592
    'ENCODING_MAX_COMPRESSION_RATIO': _LABEL,   # 178 of 3,592
    'ENCODING_MIN_COMPRESSION_RATIO': _LABEL,   # 247 of 3,592
    'RIGHT_ASCENSION'               : _LABEL,   # 1 of 14,573, 314.525 for -45.4753
    'SUB_SOLAR_LATITUDE'            : _LABEL,   # 41 of 3,880
    'SUB_SOLAR_LONGITUDE'           : _LABEL,   # 90 of 3,880
    'SUB_SPACECRAFT_LATITUDE'       : _LABEL,   # 11 of 3,753
    'SUB_SPACECRAFT_LINE'           : _LABEL,   # 134 of 3,753
    'SUB_SPACECRAFT_LINE_SAMPLE'    : _LABEL,   # 148 of 3,753
    'TARGET_CENTER_DISTANCE'        : _LABEL,   # 139 of 3,880

    # The labels correct the header: 150 respell BLACK_SKY or NON_SCIENCE with a space,
    # and 51 name a ring image's target J RINGS where the header names JUPITER; 201 of
    # 14,672.
    'TARGET_NAME'                   : _LABEL,
}

# The relative error allowed in comparing a rounding difference to half a unit
_ROUNDING_TOLERANCE = 1.e-6


def is_null(value):
    """True if a value from either source stands for no value.

    Parameters:
        value (Any): A value from the PDS3 label or the VICAR header.

    Returns:
        bool: True for a PDS3 null constant (UNK, N/A, or NULL), in any case.
    """

    return isinstance(value, str) and value.strip().upper() in _PDS3_NULLS


def label_items(label_dict):
    """The metadata items of a PDS3 label.

    Parameters:
        label_dict (dict): The label, as given by `Pds3Label.as_dict()`.

    Returns:
        dict: Every value of the label that is not a PDS3 null. A set, such as
        CUT_OUT_WINDOW, is given as the list of its values in label order. An
        ORBIT_NUMBER that is not an orbit of the Jupiter tour (0-35) is omitted.
    """

    items = {}
    for (key, value) in label_dict.items():
        if isinstance(value, set) and key + '_list' in label_dict:
            value = label_dict[key + '_list']
        if not is_null(value):
            items[key] = value

    # The labels of GO_0003-GO_0009 and GO_0016 give the ORBIT_NUMBER of every frame as
    # the first two characters of its observation ID read as a 16-bit integer, such as
    # 12613 for "E1"; 5,213 of the 14,672 products. These frames precede the Jupiter tour
    # and have no orbit number (#41).
    if items.get('ORBIT_NUMBER') not in _ORBITS:
        items.pop('ORBIT_NUMBER', None)

    return items


def header_items(header_dict):
    """The metadata items of a VICAR header, also given under their PDS3 names.

    Parameters:
        header_dict (dict): The header, as given by `VicarLabel.as_dict(unique=False)`.

    Returns:
        dict: Every value of the header that is set, with string values stripped of
        padding, plus the same values under the PDS3 keyword names, spelled as the PDS3
        label spells them. SPACECRAFT_CLOCK_START_COUNT, IMAGE_TIME,
        EARTH_RECEIVED_TIME, FILTER_NAME, FRAME_DURATION, GAIN_MODE_ID,
        EXPOSURE_TYPE, and ORBIT_NUMBER are built from the VICAR fields that encode them.
        A value the header leaves unset (a numeric -999, a PICNO of "?", a time whose
        year is 0 or less, an unknown code, or an observation ID that names no orbit of
        the Jupiter tour) gives no item.
    """

    items = {}
    for (key, value) in header_dict.items():
        if isinstance(value, str):
            value = value.strip()
        if value != _VICAR_FILL and not is_null(value):
            items[key] = value

    for (pds3_key, vicar_key) in _VICAR_PDS3_PAIRS:
        if vicar_key in items:
            items[pds3_key] = items[vicar_key]

    if items.get('IMAGE_ID') == _NO_PICNO:
        del items['IMAGE_ID']

    if 'RIM' in items and 'MOD91' in items:
        items['SPACECRAFT_CLOCK_START_COUNT'] = (f'{items["RIM"]:08d}.'
                                                 f'{items["MOD91"]:02d}')

    _set_coded(items, 'FILTER_NAME', 'FILTER', dict(enumerate(_FILTER_NAMES)))
    _set_coded(items, 'FRAME_DURATION', 'RATE', _FRAME_DURATIONS)
    _set_coded(items, 'GAIN_MODE_ID', 'GAIN', _GAIN_MODES)

    # The last flag of FIBE, or of MOFIBE in the phase-2 headers of the Jupiter tour,
    # marks an extended exposure
    flags = items.get('MOFIBE', items.get('FIBE'))
    if flags:
        items['EXPOSURE_TYPE'] = 'EXTENDED' if flags[-1] == '1' else 'NORMAL'

    for (pds3_key, prefix) in (('IMAGE_TIME', 'SCET'), ('EARTH_RECEIVED_TIME', 'ERT')):
        time = _vicar_time(items, prefix)
        if time is not None:
            items[pds3_key] = time

    orbit = _ORBIT_CODE.match(items.get('OBSERVATION_ID', ''))
    if orbit and int(orbit.group().lstrip('CEGJ')) in _ORBITS:
        items['ORBIT_NUMBER'] = int(orbit.group().lstrip('CEGJ'))

    return items


def _set_coded(items, pds3_key, vicar_key, names):
    """Give the PDS3 item that a VICAR code stands for, if the code is known.

    Parameters:
        items (dict): The header's items, modified in place.
        pds3_key (str): The name of the PDS3 item.
        vicar_key (str): The name of the VICAR keyword that holds the code.
        names (dict): The PDS3 value for each known code.
    """

    code = items.get(vicar_key)
    if code in names:
        items[pds3_key] = names[code]


def _vicar_time(items, prefix):
    """The time given by a set of VICAR time fields.

    Parameters:
        items (dict): The header's items.
        prefix (str): The prefix of the time fields, "SCET" or "ERT".

    Returns:
        datetime.datetime | None: The time, or None if the header gives none, which it
        shows by a year of 0 or less, or by leaving out the fields.
    """

    fields = [items.get(prefix + name) for name in ('YEAR', 'DAY', 'HOUR', 'MIN', 'SEC',
                                                    'MSEC')]
    if None in fields or fields[0] <= 0:
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
        has the value chosen by :func:`resolve`.
    """

    merged = dict(header or {})
    merged.update(label or {})
    if label and header:
        for key in label.keys() & header.keys():
            merged[key] = resolve(key, label[key], header[key])
    return merged


def resolve(key, label_value, header_value):
    """The value of an item that both sources give.

    A number that is the other source's number rounded loses to it. Otherwise, where the
    values differ, the item's rule names the source whose value is kept, and the PDS3
    label's value is kept for any item without a rule.

    Parameters:
        key (str): The item's name.
        label_value (Any): The PDS3 label's value.
        header_value (Any): The VICAR header's value.

    Returns:
        Any: The value kept.
    """

    if label_value == header_value:
        return label_value
    if _is_rounded(label_value, header_value):
        return header_value
    if _is_rounded(header_value, label_value):
        return label_value
    return header_value if _RULES.get(key) == _HEADER else label_value


def _is_rounded(rounded, precise):
    """True if one number is another rounded to fewer digits.

    Parameters:
        rounded (Any): The value that may be rounded.
        precise (Any): The value that may be more precise.

    Returns:
        bool: True if both are floats, `rounded` has fewer significant digits than
        `precise`, and the two differ by no more than half a unit in the last digit of
        `rounded`.
    """

    if not (isinstance(rounded, float) and isinstance(precise, float)):
        return False

    (digits, place) = _precision(rounded)
    if digits >= _precision(precise)[0]:
        return False

    # A value written to its last place in one source and rounded half-up in the other can
    # differ by exactly half a unit, to within the error of its binary representation
    return abs(precise - rounded) <= 0.5 * 10.**place * (1. + _ROUNDING_TOLERANCE)


def _precision(number):
    """The significant digits of a float as Python writes it, and its last digit's place.

    Parameters:
        number (float): The number.

    Returns:
        tuple[int, int]: The count of significant digits, which is zero for zero, and the
        power of ten of the last significant digit, which is zero for zero. The trailing
        zeros of a whole number are not counted as significant.
    """

    (mantissa, _, exponent) = repr(number).lstrip('-').partition('e')
    (whole, _, fraction) = mantissa.partition('.')
    fraction = fraction.rstrip('0')
    significant = (whole + fraction).lstrip('0')
    place = int(exponent or 0) - len(fraction)
    if not fraction:
        stripped = significant.rstrip('0')
        place += len(significant) - len(stripped)
        significant = stripped
    if not significant:
        place = 0
    return (len(significant), place)

##########################################################################################
