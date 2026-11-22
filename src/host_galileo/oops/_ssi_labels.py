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

from oops.host import is_pds3_null

# The value of a numeric VICAR keyword that the header leaves unset, such as SMRAZ; the
# indexes also write it, as their NULL_CONSTANT, for a missing number
_VICAR_FILL = -999.

# The value of a VICAR keyword that the header gives for no value. PICNO "?" is in 162
# and SPICE_C_ID "?" in 64 of the 14,672 products. The others were found in a sample of
# 880 headers, 40 from each volume: TCA "CLOSEST TIME" in 287 of the sample's 600 phase-1
# headers, where the others give a time from closest approach such as "-000T02:21:45Z";
# and SMEAR 0.1 in all 280 headers of the Jupiter tour, where every earlier header gives
# the -999 fill and every PDS3 label gives UNK for SMEAR_MAGNITUDE. The rest are the
# values that each volume's CATALOG/DATASET.CAT documents as "not applicable" (RAD, the
# ring radius, of 0.0) or "not implemented" (CONE and SMRAZ of 0.0, and SOLRANGE of
# 7.779091e8 km, which is Jupiter's distance from the Sun even in the images of Venus,
# the Earth, and Gaspra that give it). DATASET.CAT also lists, for targets other than a
# planetary body or asteroid (its section C) and for images without SPICE data (section
# D), an HSCL, VSCL, SLRANGE, and PLRANGE of 0.0, which are never real values: 3,190 of
# the 14,672 headers give HSCL and VSCL 0, 2,942 SLRANGE 0, and 2,672 PLRANGE 0.
_VICAR_FILLS = {'PICNO': '?', 'SPICE_C_ID': '?', 'TCA': 'CLOSEST TIME', 'SMEAR': 0.1,
                'RAD': 0., 'CONE': 0., 'SMRAZ': 0., 'SOLRANGE': 777909100.,
                'HSCL': 0., 'VSCL': 0., 'SLRANGE': 0., 'PLRANGE': 0.}

# The PDS3 items of the documented values that are never real, in the PDS3 labels, the
# VICAR header after its conversion to PDS3 names, and the indexes. The DATASET.CAT of
# each volume lists CENTER_RING_RADIUS = 0.0 as not applicable; CONE_ANGLE = 0.000,
# SMEAR_AZIMUTH = 0.000, and SOLAR_DISTANCE = 7.779091e8 as not implemented; and a
# CENTRAL_BODY_DISTANCE, HORIZONTAL_PIXEL_SCALE, VERTICAL_PIXEL_SCALE, and SLANT_DISTANCE
# of 0 as not available in sections C and D (3,190 of the 14,672 labels give both pixel
# scales 0, 2,942 SLANT_DISTANCE 0, and 748 CENTRAL_BODY_DISTANCE 0). The product labels
# and headers are part of the holdings and keep these values.
#
# The indexes write -999 for most of them, but 5,354 of the 10,669 index rows before the
# Jupiter tour give the SOLAR_DISTANCE 777909100 (rms-data-projects #258, item 12), and
# 248 give both pixel scales 0 where the column's null constant is -999 (rms-data-projects
# #258, item 14). Until the index is corrected, these values are dropped from index rows
# too; a corrected row writes the null constant, which is dropped as well, so the rule
# gives the same result on a corrected table.
_PDS3_FILLS = {'CENTER_RING_RADIUS': 0., 'CONE_ANGLE': 0., 'SMEAR_AZIMUTH': 0.,
               'SOLAR_DISTANCE': 777909100., 'CENTRAL_BODY_DISTANCE': 0.,
               'HORIZONTAL_PIXEL_SCALE': 0., 'VERTICAL_PIXEL_SCALE': 0.,
               'SLANT_DISTANCE': 0.}

# The items, under their PDS3 and VICAR names, that are never negative. A survey of all
# 14,672 labels and headers found 128 GO_0017 headers with HSCL = VSCL = -999000, two of
# whose labels (C3/RINGS/C0368976000R and C0368992339R) give HORIZONTAL_PIXEL_SCALE =
# VERTICAL_PIXEL_SCALE = -9.99e5, and one GO_0018 product (E4/EUROPA/C0374667300R) whose
# label gives SLANT_DISTANCE = -547825 and header SLRANGE = -547825, none of them
# declared as a fill (rms-data-projects #259, item 11). A negative value of these items is
# no value. The indexes write these columns as E12.6 and cut the last digit of a negative
# value with a two-digit exponent, so their rows for the same products give both pixel
# scales -9.99 (2 rows of GO_0017, and the same 2 of GO_0999) and SLANT_DISTANCE -5.47825
# (1 row of GO_0018, and the same row of GO_0999) (rms-data-projects #258, item 15).
# These are the only negative values with a cut exponent in the 23 indexes. Until the
# index is
# corrected, the same rule drops those values; a corrected row writes the null constant
# -999, which is dropped as well.
_NEVER_NEGATIVE = ('HORIZONTAL_PIXEL_SCALE', 'VERTICAL_PIXEL_SCALE', 'SLANT_DISTANCE',
                   'HSCL', 'VSCL', 'SLRANGE')

# The groups of items that DATASET.CAT lists as 0 where the value is not available: the
# pointing, when the C kernels had none for the image (section D; 93 labels and headers);
# the latitude and longitude of the image center (CENTER_LATITUDE and CENTER_LONGITUDE in
# a label, LAT and LON in a header, which the host gives as INTERCEPT_POINT_LATITUDE and
# INTERCEPT_POINT_LONGITUDE), both 0 in 4,028 labels and headers, including 1,030 with a
# real pixel scale whose image center misses the target; and the three azimuths, all 0 in
# the same 3,190 labels and headers whose pixel scales are 0. One item of a group can be
# 0 in a real value, so a group counts as no value only when every item of it that the
# source gives is 0 and the source gives at least two of them. The indexes, which lack
# CENTER_LATITUDE, CENTER_LONGITUDE, and SUB_SPACECRAFT_AZIMUTH, write -999 for these
# items where the labels write 0.
_ZERO_GROUPS = (
    ('RIGHT_ASCENSION', 'DECLINATION', 'TWIST_ANGLE'),
    ('CENTER_LATITUDE', 'CENTER_LONGITUDE'),
    ('INTERCEPT_POINT_LATITUDE', 'INTERCEPT_POINT_LONGITUDE'),
    ('NORTH_AZIMUTH', 'SUB_SOLAR_AZIMUTH', 'SUB_SPACECRAFT_AZIMUTH'),
)

# The times from closest approach that DATASET.CAT lists as "-73T..." where the value is
# not available (sections C and D). 48 labels give "-73T69:12:71Z", which is not a
# valid time, as TIME_FROM_CLOSEST_APPROACH and 48 as SATELLITE_TIME_FROM_CLST_APR, and
# 48 index rows give it in each of the two columns (rms-data-projects #258, item 14); a
# real time has three digits of days, such as "-023T13:24:21Z". Until the
# index is corrected, the index value is dropped too; a corrected row writes a null,
# which is dropped as well.
_UNKNOWN_TIME_ITEMS = ('NTV_TIME_FROM_CLOSEST_APPROACH', 'NTV_SAT_TIME_FROM_CLOSEST_APRH')
_UNKNOWN_TIME_PREFIX = '-73T'

# The labels before the Jupiter tour give the time from closest approach under other names
# than the tour labels, the indexes, and the VICAR header's TCA do; the host gives them
# under the tour labels' names
_LABEL_NAMES = {
    'TIME_FROM_CLOSEST_APPROACH'  : 'NTV_TIME_FROM_CLOSEST_APPROACH',
    'SATELLITE_TIME_FROM_CLST_APR': 'NTV_SAT_TIME_FROM_CLOSEST_APRH',
}

# A value that an index or supplemental index column declares as its null or
# not-applicable constant is no value, and the table reader marks it in the row's "_mask"
# key. The SSI tables declare the supplement's CUT_OUT_WINDOW and TRUTH_WINDOW of -1, -1,
# -1, -1, its ICT_DESPIKE_THRESHOLD and ICT_QUANTIZATION_STEP_SIZE of -1, the index's
# compression ratios of 0, its numeric -999, and its ORBIT_NUMBER of -99. Comparing the
# 14,662 rows of GO_0002-GO_0023 with both labels of each product: in 11,080 rows that
# hold one of the window, ICT, or ratio constants, neither label gives the item at all; in
# the other 646 (TRUTH_WINDOW and the ICT values) or 647-662 (the compression ratios), the
# labels give a real value that the index replaced, such as an ICT_DESPIKE_THRESHOLD of
# 255 or a compression ratio of 1.831, and only 1-3 rows per ratio hold 0 where the labels
# also write 0. Until the index is corrected (rms-data-projects #258, item 1), those
# values are lost. A corrected table declares the constants that it uses, so the host
# follows the declarations rather than these values.
#
# The tables also declare three real values, which the product labels give too, as the
# column's null or not-applicable constant: EXPOSURE_DURATION 0 (301 rows; NULL_CONSTANT),
# ON_CHIP_MOSAIC_FLAG N, and PRODUCT_VERSION_ID 0 (NOT_APPLICABLE_CONSTANT). Until the
# index is corrected (rms-data-projects #258, item 2), the host keeps each of these
# values in its column even though the table reader masks it. Only that value is kept:
# a corrected table that declares another constant in one of these columns, such as a
# PRODUCT_VERSION_ID of -1, has that constant dropped, and one that no longer declares
# these values gives them unmasked.
_REAL_NULL_CONSTANTS = {'EXPOSURE_DURATION': 0., 'ON_CHIP_MOSAIC_FLAG': 'N',
                        'PRODUCT_VERSION_ID': 0}

# The value that an index column writes where it has none, without declaring it. The
# index's SMEAR_MAGNITUDE is 0.1 in all 3,993 rows of the Jupiter tour and -999 in every
# earlier row, where every label gives UNK; the 0.1 is the VICAR header's placeholder
# (_VICAR_FILLS) and never a real value (rms-data-projects #258, item 11).
_INDEX_PLACEHOLDERS = {'SMEAR_MAGNITUDE': 0.1}

# The index columns that hold the three compression ratios in a table that has them in
# rotated columns, and the ratio each holds. Every current SSI index label names the
# three columns, in the order of their bytes in the row, ENCODING_MIN_COMPRESSION_RATIO,
# ENCODING_MAX_COMPRESSION_RATIO, and ENCODING_COMPRESSION_RATIO, but the values in them
# are each image's ratio, minimum, and maximum. So the index gives each image's maximum
# under ENCODING_COMPRESSION_RATIO, its ratio under the minimum's name, and its minimum
# under the maximum's name: in all 2,935 rows with a nonzero ratio (GO_0017-GO_0023 and
# the cumulative GO_0999), the three values are the product label's in that rotated
# order, as in GO_0021/C21/IO/C0506406074R (index 5.428, 3.758, 3.03; label 3.758, 3.03,
# 5.428). Until the index is corrected (rms-data-projects #258, item 8), the host gives
# each value under the name of the ratio it holds, in a table whose label shows the
# defect (ratios_rotated).
_INDEX_RATIO_COLUMNS = {
    'ENCODING_COMPRESSION_RATIO'    : 'ENCODING_MAX_COMPRESSION_RATIO',
    'ENCODING_MIN_COMPRESSION_RATIO': 'ENCODING_COMPRESSION_RATIO',
    'ENCODING_MAX_COMPRESSION_RATIO': 'ENCODING_MIN_COMPRESSION_RATIO',
}

# The ratio columns in the order of their bytes in a table whose label shows the defect.
# The decision rests on this order alone. The correction asked for in rms-data-projects
# #258, item 8, renames the three columns and leaves the values in place, which changes
# this order.
_ROTATED_RATIO_ORDER = ('ENCODING_MIN_COMPRESSION_RATIO',
                        'ENCODING_MAX_COMPRESSION_RATIO',
                        'ENCODING_COMPRESSION_RATIO')

# The index column names that differ from the product labels' names for the same item:
# the indexes drop the last letter of the tour labels' NTV_SAT_TIME_FROM_CLOSEST_APRH.
# Until the index is corrected (rms-data-projects #258, item 13), the host gives the
# labels' name; a corrected table already uses it.
_INDEX_NAMES = {'NTV_SAT_TIME_FROM_CLOSEST_APR': 'NTV_SAT_TIME_FROM_CLOSEST_APRH'}

# The index values that the product labels give differently, and the labels' value. The
# indexes give all 3,993 images of the Jupiter tour data set version 1.0, where every
# product label gives version 1.1. Until the index is
# corrected (rms-data-projects #258, item 9), the host gives the labels' version; a
# corrected row already gives it.
_INDEX_VALUE_REPAIRS = {
    'DATA_SET_ID': {'GO-J/JSA-SSI-2-REDR-V1.0': 'GO-J/JSA-SSI-2-REDR-V1.1',
                    'GO-J/JSA-SSI-4-REDR-V1.0': 'GO-J/JSA-SSI-4-REDR-V1.1'},
}

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
# encounter target followed by its single digit. The IDs of the earlier mission phases
# begin with a code of that phase instead, such as "E1" and "E2" for the Earth flybys,
# "GA" for Gaspra, and "VE" for Venus.
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


def label_json(label_dict):
    """The PDS3 label's dictionary in a form that `json` can write.

    Parameters:
        label_dict (dict): The label, as given by `Pds3Label.as_dict()`, or one of its
            objects.

    Returns:
        dict: Every keyword of `label_dict` in the same order and with the same value,
        except that a date or time is given as the ISO 8601 text of the label's value,
        with the label's digits but without a trailing "Z" (its "_fmt" key, such as
        "1999-10-11T04:29:52.510"), or else in ISO 8601 form, a set as its "_list" key's
        list in label order or else sorted by each element's text, and an object as this
        same dictionary. A keyword whose value is None, such as END, is left out, and so
        is every key that the parser adds beside the keywords, which are those with a
        lower-case suffix ("_fmt", "_list", "_offset", "_unit", "_day", "_sec",
        "_unwrap", "_source", ...) and "objects". Of these, only a pointer's "_offset"
        holds anything that the keyword's own value does not: the record of the image
        file at which the object starts, which the host never uses (the array's place
        comes from the VICAR header), and which 52 GO_0019 labels give wrongly.
    """

    json_dict = {}
    for (key, value) in label_dict.items():
        if value is None or key != key.upper():
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
        dict: Every keyword of the label whose value is not a PDS3 null (as
        `oops.host.is_pds3_null()` decides), including each object with the same
        treatment. A null element of a list is left out: the lists of an SSI label whose
        positions carry meaning (CUT_OUT_WINDOW, TRUTH_WINDOW, and STAR_WINDOW) never
        hold one. TIME_FROM_CLOSEST_APPROACH and SATELLITE_TIME_FROM_CLST_APR are given
        under the tour labels' names, NTV_TIME_FROM_CLOSEST_APPROACH and
        NTV_SAT_TIME_FROM_CLOSEST_APRH. The values that the archive documents as standing
        for no value are left out, as :func:`drop_fills` describes, and so is an
        ORBIT_NUMBER that is not an orbit of the Jupiter tour (0-35).
    """

    items = {}
    for (key, value) in label_dict.items():
        if isinstance(value, dict):
            value = label_items(value)
        elif isinstance(value, list) and not is_pds3_null(value):
            value = [item for item in value if not is_pds3_null(item)]
        if not is_pds3_null(value):
            items[_LABEL_NAMES.get(key, key)] = value
    drop_fills(items)

    # The labels of GO_0003-GO_0009 and GO_0016 give the ORBIT_NUMBER of every frame as
    # the first two characters of its observation ID read as a 16-bit integer, such as
    # 12613 for "E1"; 5,213 of the 14,672 products. These frames precede the Jupiter tour
    # and have no orbit number (#41).
    if items.get('ORBIT_NUMBER') not in _ORBITS:
        items.pop('ORBIT_NUMBER', None)

    return items


def drop_fills(items):
    """Remove the values that the archive documents as standing for no value, and the
    negative pixel scales and slant distances.

    The same rules apply to the items of a PDS3 label, of a VICAR header once its values
    have their PDS3 names, and of an index row.

    Parameters:
        items (dict): The items of one source under their PDS3 names (and a VICAR
            header's other items under their own names), modified in place. These are
            removed: a CONE_ANGLE, SMEAR_AZIMUTH, CENTER_RING_RADIUS,
            CENTRAL_BODY_DISTANCE, HORIZONTAL_PIXEL_SCALE, VERTICAL_PIXEL_SCALE, or
            SLANT_DISTANCE of 0; a negative HORIZONTAL_PIXEL_SCALE, VERTICAL_PIXEL_SCALE,
            SLANT_DISTANCE, HSCL, VSCL, or SLRANGE; a SOLAR_DISTANCE of 777909100; an
            NTV_TIME_FROM_CLOSEST_APPROACH or NTV_SAT_TIME_FROM_CLOSEST_APRH that begins
            with "-73T"; and each of these groups when the items contain at least two of
            its members and every member they contain is 0: RIGHT_ASCENSION,
            DECLINATION, and TWIST_ANGLE; CENTER_LATITUDE and CENTER_LONGITUDE;
            INTERCEPT_POINT_LATITUDE and INTERCEPT_POINT_LONGITUDE; NORTH_AZIMUTH,
            SUB_SOLAR_AZIMUTH, and SUB_SPACECRAFT_AZIMUTH. A single 0 of a group is kept.
    """

    for key in [key for (key, value) in items.items() if _is_fill(key, value)]:
        del items[key]

    for group in _ZERO_GROUPS:
        given = [key for key in group if key in items]
        if len(given) >= 2 and all(items[key] == 0 for key in given):
            for key in given:
                del items[key]


def _is_fill(key, value):
    """True if an item's value is one that the archive documents as no value, or a
    negative value of an item that is never negative.

    Parameters:
        key (str): The item's PDS3 name, or a VICAR header item's own name.
        value (Any): Its value.

    Returns:
        bool: True for a value of `_PDS3_FILLS`, for a negative number in an item of
        `_NEVER_NEGATIVE`, and for a time from closest approach that begins with "-73T".
    """

    if key in _UNKNOWN_TIME_ITEMS:
        return isinstance(value, str) and value.startswith(_UNKNOWN_TIME_PREFIX)
    if key in _NEVER_NEGATIVE and isinstance(value, (int, float)) and value < 0:
        return True
    return key in _PDS3_FILLS and value == _PDS3_FILLS[key]


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
        observation ID. A value the header leaves unset is left out: a numeric -999, a
        PICNO or SPICE_C_ID of "?", a TCA of "CLOSEST TIME", a SMEAR of 0.1, the values
        that the archive documents as not implemented, not applicable, or not available
        (a RAD, CONE, SMRAZ, HSCL, VSCL, SLRANGE, or PLRANGE of 0 and a SOLRANGE of
        777909100, and, under their PDS3 names, the values that :func:`drop_fills`
        removes, such as RA, DEC, and TWIST when all three are 0, and LAT and LON when
        both are 0), a negative HSCL, VSCL, or SLRANGE, a placeholder text such as NOT
        APPLICABLE, and a set of time fields
        that is not a valid time. ORBIT_NUMBER is given only for an observation ID that
        begins with an orbit of the Jupiter tour.
    """

    items = {}
    for (key, value) in header_dict.items():
        if isinstance(value, str):
            value = value.strip()
        if (value != _VICAR_FILL and value != _VICAR_FILLS.get(key)
                and not is_pds3_null(value)):
            items[key] = value

    for (pds3_key, vicar_key) in _VICAR_PDS3_PAIRS:
        if vicar_key in items:
            items[pds3_key] = items.pop(vicar_key)
    drop_fills(items)

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


def ratios_rotated(label):
    """True if an index table writes the three compression ratios in rotated columns.

    The decision is made once for the whole table, from the order of the column names in
    its label alone: the table has the defect if its label names the three ratio
    columns, in the order of their bytes in the row, ENCODING_MIN_COMPRESSION_RATIO,
    ENCODING_MAX_COMPRESSION_RATIO, and ENCODING_COMPRESSION_RATIO, as every SSI index
    label does. A table whose label names the
    columns in another order is taken to be corrected, whenever the label was made.

    Parameters:
        label (Pds3Label | dict): The label of the index.

    Returns:
        bool: True if the table's ratio columns are rotated; False if they are in
        another order or if the label lacks one of the three columns.
    """

    for table in label.values():
        if isinstance(table, dict) and all(isinstance(table.get(name), dict)
                                           for name in _ROTATED_RATIO_ORDER):
            starts = [table[name].get('START_BYTE', 0) for name in _ROTATED_RATIO_ORDER]
            return starts == sorted(starts)
    return False


def row_items(row_dict, *, rotated_ratios=False):
    """The metadata items of a row of an SSI index and its supplement.

    Parameters:
        row_dict (dict): The row, as given by `oops.Host._index_rows()`.
        rotated_ratios (bool, optional): True if the row comes from a table that writes
            the three compression ratios in rotated columns, as :func:`ratios_rotated`
            decides from the table's label.

    Returns:
        dict: Every column value of the row as a Python str (stripped of padding), int,
        float, bool, or list. A value that stands for no value is left out: a value that
        the table reader masks because the table declares it as the column's null or
        not-applicable constant (such as the CUT_OUT_WINDOW of -1, -1, -1, -1 and the
        compression ratio of 0), except an EXPOSURE_DURATION of 0, an
        ON_CHIP_MOSAIC_FLAG of N, and a PRODUCT_VERSION_ID of 0, which are real values;
        a PDS3 null; a numeric -999 (the indexes' null constant for a missing number);
        the SMEAR_MAGNITUDE of 0.1; the values that the archive documents as standing for
        no value, as :func:`drop_fills` describes (such as the SOLAR_DISTANCE of
        777909100 and a pixel scale of 0); and an ORBIT_NUMBER that is not an orbit of
        the Jupiter tour (0-35), which the index writes as -99 for the earlier frames.
        With `rotated_ratios`, the three compression ratios are given under the names of
        the ratios they hold: the value of the index's ENCODING_COMPRESSION_RATIO as
        ENCODING_MAX_COMPRESSION_RATIO, of its ENCODING_MIN_COMPRESSION_RATIO as
        ENCODING_COMPRESSION_RATIO, and of its ENCODING_MAX_COMPRESSION_RATIO as
        ENCODING_MIN_COMPRESSION_RATIO. A DATA_SET_ID of the Jupiter tour with version
        1.0 is given as version 1.1, as every product label gives it, and
        NTV_SAT_TIME_FROM_CLOSEST_APR is given under the labels' name,
        NTV_SAT_TIME_FROM_CLOSEST_APRH. The "_mask" keys themselves are left out. A row
        that lacks them, as rows built by hand may, keeps the values that they would
        have masked.
    """

    items = {}
    for (column, value) in row_dict.items():
        if column.endswith('_mask'):
            continue
        value = value.tolist() if isinstance(value, (np.ndarray, np.generic)) else value
        if isinstance(value, str):
            value = value.strip()
        if (np.all(row_dict.get(column + '_mask', False))
                and not (column in _REAL_NULL_CONSTANTS
                         and value == _REAL_NULL_CONSTANTS[column])):
            continue
        key = _INDEX_RATIO_COLUMNS.get(column, column) if rotated_ratios else column
        key = _INDEX_NAMES.get(key, key)
        if isinstance(value, str):
            value = _INDEX_VALUE_REPAIRS.get(key, {}).get(value, value)
        if (value != _VICAR_FILL and value != _INDEX_PLACEHOLDERS.get(key)
                and not is_pds3_null(value)):
            items[key] = value
    drop_fills(items)

    if 'ORBIT_NUMBER' in items and items['ORBIT_NUMBER'] not in _ORBITS:
        del items['ORBIT_NUMBER']

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
