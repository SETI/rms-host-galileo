##########################################################################################
# host_galileo/oops/_ssi_labels.py
##########################################################################################
"""The metadata of an SSI image, from its PDS3 label, its VICAR header, or both.

Each source is first reduced to its metadata items. The PDS3 label gives its values under
their own names. The VICAR header gives each value that the host converts into a PDS3
item only under the PDS3 name, spelled as the label spells it, and every other value under
its own name. A value that a source gives as "no value" (a PDS3 null such as UNK, or a
VICAR fill such as PICNO "?" or -999) is not an item of that source. Where both sources
give an item, :func:`merge_items` keeps the value of the item's one fixed source.

Every value is in a form that `json` can write: times are ISO 8601 strings as the PDS3
label writes them, and sets are lists.
"""

import calendar
import datetime as dt
import re

import numpy as np

# The PDS3 constants that stand for a missing value
_PDS3_NULLS = ('UNK', 'N/A', 'NULL')

# The value of a numeric VICAR keyword that the header leaves unset, such as SMRAZ; the
# indexes also write it, as their NULL_CONSTANT, for a missing number
_VICAR_FILL = -999.

# The value of a VICAR string keyword that the header leaves unset: PICNO in 162 and
# SPICE_C_ID in 64 of the 14,672 products
_VICAR_STRING_FILLS = {'PICNO': '?', 'SPICE_C_ID': '?'}

# The CUT_OUT_WINDOW value written in the supplemental index when the label has none
_NO_WINDOW = [-1, -1, -1, -1]

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

# The valid range of each VICAR time field; the day of the year is also checked against
# the year's length. A year of 0 or less is the header's fill (-32768 in the RAW_CAL
# frames), and a second of 60 is a leap second.
_TIME_LIMITS = {'YEAR': (1, dt.MAXYEAR), 'DAY': (1, 366), 'HOUR': (0, 23),
                'MIN': (0, 59), 'SEC': (0, 60), 'MSEC': (0, 999)}

# The orbits of the Jupiter tour, as the code that begins an observation ID: the orbit
# number as two digits from orbit 10 on, and before that the initial of the orbit's
# encounter target followed by its single digit. The IDs of the earlier mission phases,
# including those of the Earth flybys, which begin "E1" and "E2", name no orbit.
_ORBIT_CODE = re.compile(r'\d\d|J0|G1|G2|C3|E4|J5|E6|G7|G8|C9')

# The orbit numbers of the Jupiter tour, from J0 to orbit 35
_ORBITS = range(36)

# The resolution rules: for each item that the two sources give differently, the one
# source whose value is kept for every image. An item without a rule keeps the PDS3
# label's value. A value missing from the winning source is not a difference; the other
# source's value is used. The counts come from comparing the two labels of all 14,672 SSI
# products of GO_0002-GO_0023 ("of" counts the products that give the item in both), and
# the SPICE checks evaluate the host's own kernels at each image's IMAGE_TIME. Every rule,
# with its evidence and the images on which its source is wrong, is recorded in #41.
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
    'NORTH_AZIMUTH'                 : _HEADER,  # 6,948 of 14,550
    'PHASE_ANGLE'                   : _HEADER,  # 8,720 of 10,545
    'SUB_SOLAR_AZIMUTH'             : _HEADER,  # 213 of 14,550
    'SUB_SPACECRAFT_AZIMUTH'        : _HEADER,  # 4,898 of 14,550
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
    # of 1,748 tour products, where the label is off by a median 64 km. 8,745 of 12,648
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


def is_null(value):
    """True if a value from either source, or from an index row, stands for no value.

    Parameters:
        value (Any): A value from the PDS3 label, the VICAR header, or an index.

    Returns:
        bool: True for None, for a PDS3 null constant (UNK, N/A, or NULL) in any case and
        with any padding, and for a non-empty list of nothing but these.
    """

    if isinstance(value, list):
        return len(value) > 0 and all(is_null(item) for item in value)
    return value is None or (isinstance(value, str)
                             and value.strip().upper() in _PDS3_NULLS)


def label_json(label_dict):
    """The PDS3 label's dictionary in a form that `json` can write.

    Parameters:
        label_dict (dict): The label, as given by `Pds3Label.as_dict()`, or one of its
            objects.

    Returns:
        dict: Every key of `label_dict` in the same order and with the same value, except
        that a date or time is given as the ISO 8601 text of the label's value, with the
        label's digits but without a trailing "Z" (its "_fmt" key, such as
        "1999-10-11T04:29:52.510"), or else in ISO 8601 form, a set as its "_list" key's
        list in label order or else sorted by each element's text, and an object as this
        same dictionary. A key whose value is None, such as END, and the label text that
        the parser gives for a keyword (its "_source" key) are left out.
    """

    json_dict = {}
    for (key, value) in label_dict.items():
        if value is None or (key.endswith('_source') and key[:-7] in label_dict):
            continue
        if isinstance(value, (dt.date, dt.time)):
            value = label_dict.get(key + '_fmt', value.isoformat())
        elif isinstance(value, set):
            # A set can mix numbers and text, as the parser reads some STAR_WINDOW sets
            value = (label_dict[key + '_list'] if key + '_list' in label_dict
                     else sorted(value, key=str))
        elif isinstance(value, dict):
            value = label_json(value)
        json_dict[key] = value
    return json_dict


def header_json(label):
    """The VICAR header's dictionary, with every parameter, in a form `json` can write.

    Parameters:
        label (VicarLabel): The VICAR header.

    Returns:
        dict: Every parameter of the header in header order. A parameter that occurs once
        is keyed by its name; each occurrence of a repeated one, such as the TASK, USER,
        and DAT_TIM of each step of the processing history, by its name followed by its
        occurrence counted from zero in brackets, such as "TASK[1]".
    """

    return {(key if isinstance(key, str) else f'{key[0]}[{key[1]}]'): value
            for (key, value) in label.as_dict(unique=True).items()}


def label_items(label_dict):
    """The metadata items of a PDS3 label.

    Parameters:
        label_dict (dict): The label, as given by :func:`label_json`.

    Returns:
        dict: Every keyword of the label whose value is not a PDS3 null, including each
        object with the same treatment, and none of the keys that the parser adds (those
        with a lower-case suffix such as "_offset", "_fmt", or "_list", and "objects").
        A null element of a list is left out: the lists of an SSI label whose positions
        carry meaning (CUT_OUT_WINDOW, TRUTH_WINDOW, and STAR_WINDOW) never hold one. An
        ORBIT_NUMBER that is not an orbit of the Jupiter tour (0-35) is omitted.
    """

    items = {}
    for (key, value) in label_dict.items():
        if key != key.upper():
            continue
        if isinstance(value, dict):
            value = label_items(value)
        elif isinstance(value, list) and not is_null(value):
            value = [item for item in value if not is_null(item)]
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
    """The metadata items of a VICAR header, under PDS3 names where converted.

    Parameters:
        header_dict (dict): The header, as given by :func:`header_json`.

    Returns:
        dict: Every value of the header that is set, with string values stripped of
        padding. A value that the host converts into a PDS3 item is given only under the
        PDS3 name, spelled as the PDS3 label spells it: the keywords that name the same
        value (such as EXP for EXPOSURE_DURATION, FILTER for FILTER_NUMBER, and PA for
        OBSERVATION_ID), RIM and MOD91 as SPACECRAFT_CLOCK_START_COUNT, the SCET and ERT
        fields as IMAGE_TIME and EARTH_RECEIVED_TIME (each an ISO 8601 string such as
        "1996-06-26T08:46:04.890"), and the codes RATE and GAIN as FRAME_DURATION and
        GAIN_MODE_ID. FILTER_NAME comes from the filter number. MOD10 and MOD8, and the
        flags FIBE and MOFIBE, which give EXPOSURE_TYPE but more besides, keep their own
        names, as does a code with no PDS3 value. ORBIT_NUMBER comes from the
        observation ID. A value the header leaves unset (a numeric -999, a PICNO or
        SPICE_C_ID of "?", a set of time fields that is not a valid time, or an
        observation ID that names no orbit of the Jupiter tour) gives no item.
    """

    items = {}
    for (key, value) in header_dict.items():
        if isinstance(value, str):
            value = value.strip()
        if (value != _VICAR_FILL and value != _VICAR_STRING_FILLS.get(key)
                and not is_null(value)):
            items[key] = value

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

    orbit = _ORBIT_CODE.match(items.get('OBSERVATION_ID', ''))
    if orbit and int(orbit.group().lstrip('CEGJ')) in _ORBITS:
        items['ORBIT_NUMBER'] = int(orbit.group().lstrip('CEGJ'))

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
    """The time given by a set of VICAR time fields, as a PDS3 label writes it.

    Parameters:
        fields (list[int | None]): The year, day of the year, hour, minute, second, and
            millisecond, each None if the header leaves it out.

    Returns:
        str | None: The time as "YYYY-MM-DDTHH:MM:SS.sss", or None if the header gives
        none, which it shows by leaving out a field, by a year of 0 or less, or by a
        field out of its range, such as an hour of 24 or a day beyond the year's length.
    """

    if None in fields:
        return None
    for (value, name) in zip(fields, _TIME_FIELDS, strict=True):
        (low, high) = _TIME_LIMITS[name]
        if not low <= value <= high:
            return None

    (year, doy, hour, minute, second, msec) = fields
    if doy > (366 if calendar.isleap(year) else 365):
        return None
    date = dt.date(year, 1, 1) + dt.timedelta(days=doy - 1)
    return f'{date.isoformat()}T{hour:02d}:{minute:02d}:{second:02d}.{msec:03d}'


def row_items(row_dict):
    """The metadata items of a row of an SSI index and its supplement.

    Parameters:
        row_dict (dict): The row, as given by `oops.Host._read_index_rows()`.

    Returns:
        dict: Every column value of the row as a Python str (stripped of padding), int,
        float, bool, or list, except a PDS3 null, a numeric -999 (the indexes' null
        constant for a missing number), and the CUT_OUT_WINDOW of -1, -1, -1, -1 that the
        supplement writes for an image with no cutout window. The "_mask" keys that the
        table reader adds are left out, and its masks are not applied: the SSI indexes
        declare real values as null constants, such as an EXPOSURE_DURATION of 0 and an
        ON_CHIP_MOSAIC_FLAG of N.
    """

    items = {}
    for (key, value) in row_dict.items():
        if key.endswith('_mask'):
            continue
        value = value.tolist() if isinstance(value, (np.ndarray, np.generic)) else value
        if isinstance(value, str):
            value = value.strip()
        if value != _VICAR_FILL and not is_null(value):
            items[key] = value

    if items.get('CUT_OUT_WINDOW') == _NO_WINDOW:
        del items['CUT_OUT_WINDOW']

    return items


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

    Where the values differ, the item's rule names the one source whose value is kept for
    every image; the PDS3 label's value is kept for an item without a rule.

    Parameters:
        key (str): The item's name.
        label_value (Any): The PDS3 label's value.
        header_value (Any): The VICAR header's value.

    Returns:
        Any: The value kept.
    """

    return header_value if _RULES.get(key) == _HEADER else label_value

##########################################################################################
