##########################################################################################
# tests/oops/test_ssi_labels.py
##########################################################################################
"""Tests for the merging of an SSI image's two labels, host_galileo.oops._ssi_labels."""

import datetime as dt
import json
from typing import Any

import numpy as np
import pytest
import vicar

from host_galileo.oops import _ssi_labels


##########################################################################################
# label_json() and header_json()
##########################################################################################


def test_label_json_writes_a_time_as_the_label_does() -> None:
    raw = {'IMAGE_TIME': dt.datetime(1999, 10, 11, 4, 29, 52, 510000),
           'IMAGE_TIME_fmt': '1999-10-11T04:29:52.510'}
    assert _ssi_labels.label_json(raw)['IMAGE_TIME'] == '1999-10-11T04:29:52.510'


def test_label_json_writes_a_time_without_its_text_in_iso_form() -> None:
    raw = {'START': dt.date(1999, 10, 11)}
    assert _ssi_labels.label_json(raw)['START'] == '1999-10-11'


def test_label_json_writes_a_set_in_label_order() -> None:
    raw = {'CUT_OUT_WINDOW': {1, 400, 800}, 'CUT_OUT_WINDOW_list': [1, 1, 400, 800]}
    assert _ssi_labels.label_json(raw)['CUT_OUT_WINDOW'] == [1, 1, 400, 800]


def test_label_json_writes_a_set_without_its_list_sorted() -> None:
    assert _ssi_labels.label_json({'X': {'B', 'A'}})['X'] == ['A', 'B']


# The parser reads some STAR_WINDOW sets of GO_0017 C3/OPNAV as a mix of numbers and text
def test_label_json_writes_a_mixed_set() -> None:
    raw = {'X': {801, 'N/A'}, 'X_list': [801, 'N/A'], 'Y': {2, 'A'}}
    assert _ssi_labels.label_json(raw) == {'X': [801, 'N/A'], 'Y': [2, 'A']}


def test_label_json_drops_none() -> None:
    assert _ssi_labels.label_json({'A': 1, 'END': None}) == {'A': 1}


def test_label_json_drops_the_parsers_keys() -> None:
    raw = {'^IMAGE': 'X.IMG', '^IMAGE_offset': 7, '^IMAGE_unit': '',
           '^IMAGE_fmt': '("X.IMG", 7)', 'IMAGE_TIME': dt.datetime(2000, 1, 1),
           'IMAGE_TIME_day': 0, 'IMAGE_TIME_sec': 43200., 'IMAGE_TIME_fmt': '2000-01-01',
           'NOTE_unwrap': 'A B', 'A_source': 'A = 1', 'DATA_SOURCE': 'X',
           'objects': ['IMAGE']}
    assert _ssi_labels.label_json(raw) == {'^IMAGE': 'X.IMG', 'IMAGE_TIME': '2000-01-01',
                                           'DATA_SOURCE': 'X'}


def test_label_json_converts_an_object() -> None:
    raw = {'IMAGE': {'OBJECT': 'IMAGE', 'T': dt.time(1, 2, 3), 'END': None}}
    assert _ssi_labels.label_json(raw) == {'IMAGE': {'OBJECT': 'IMAGE', 'T': '01:02:03'}}


def test_label_json_is_serializable() -> None:
    raw = {'T': dt.datetime(2000, 1, 1), 'S': {1, 2}, 'P': ('X.IMG', 3), 'END': None}
    json.dumps(_ssi_labels.label_json(raw))


def test_header_json_keeps_every_repeated_parameter() -> None:
    label = vicar.VicarLabel([('LBLSIZE', 0), ('TASK', 'A'), ('USER', 'U'),
                              ('TASK', 'B'), ('MISSION', 'GALILEO')])
    header = _ssi_labels.header_json(label)
    assert {key: header[key] for key in ('TASK[0]', 'TASK[1]', 'MISSION')} == {
        'TASK[0]': 'A', 'TASK[1]': 'B', 'MISSION': 'GALILEO'}


def test_header_json_names_a_single_parameter_plainly() -> None:
    label = vicar.VicarLabel([('LBLSIZE', 0), ('TASK', 'A')])
    assert 'TASK' in _ssi_labels.header_json(label)

##########################################################################################
# label_items()
##########################################################################################


def test_label_items_drop_nulls() -> None:
    assert _ssi_labels.label_items({'SMEAR_AZIMUTH': 'UNK', 'X': 1}) == {'X': 1}


def test_label_items_drop_null_list_elements() -> None:
    label = {'SOURCE_PRODUCT_ID': ['A.BSP', 'N/A', 'B.BC', 'NULL']}
    assert _ssi_labels.label_items(label)['SOURCE_PRODUCT_ID'] == ['A.BSP', 'B.BC']


def test_label_items_drop_a_list_of_nulls() -> None:
    assert _ssi_labels.label_items({'SPICE_FILE_NAME': ['UNK']}) == {}


@pytest.mark.parametrize(('orbit', 'kept'), [(0, True), (35, True), (36, False),
                                             (12613, False), ('E1', False)])
def test_label_items_keep_only_tour_orbits(orbit: Any, kept: bool) -> None:
    assert ('ORBIT_NUMBER' in _ssi_labels.label_items({'ORBIT_NUMBER': orbit})) is kept


# The values that DATASET.CAT documents as not implemented, not applicable, or not
# available
_PDS3_FILLS = [('CONE_ANGLE', 0.), ('SMEAR_AZIMUTH', 0.), ('CENTER_RING_RADIUS', 0.),
               ('SOLAR_DISTANCE', 777909100.), ('CENTRAL_BODY_DISTANCE', 0.),
               ('HORIZONTAL_PIXEL_SCALE', 0.), ('VERTICAL_PIXEL_SCALE', 0.),
               ('SLANT_DISTANCE', 0.)]


@pytest.mark.parametrize(('key', 'fill'), _PDS3_FILLS)
def test_label_items_documented_fill_is_no_value(key: str, fill: float) -> None:
    assert _ssi_labels.label_items({key: fill}) == {}


@pytest.mark.parametrize('label', [{'SOLAR_DISTANCE': 7.41e8},
                                   {'HORIZONTAL_PIXEL_SCALE': 2024.28}],
                         ids=['solar-distance', 'pixel-scale'])
def test_label_items_keep_a_real_value(label: dict[str, Any]) -> None:
    assert _ssi_labels.label_items(label) == label


# The groups of items that DATASET.CAT documents as all 0 where the value is not available
_ZERO_GROUPS = [
    {'RIGHT_ASCENSION': 0., 'DECLINATION': 0., 'TWIST_ANGLE': 0.},
    {'CENTER_LATITUDE': 0., 'CENTER_LONGITUDE': 0.},
    {'INTERCEPT_POINT_LATITUDE': 0., 'INTERCEPT_POINT_LONGITUDE': 0.},
    {'NORTH_AZIMUTH': 0., 'SUB_SOLAR_AZIMUTH': 0., 'SUB_SPACECRAFT_AZIMUTH': 0.},
]
_ZERO_GROUP_IDS = ['pointing', 'center', 'intercept', 'azimuths']


@pytest.mark.parametrize('group', _ZERO_GROUPS, ids=_ZERO_GROUP_IDS)
def test_label_items_group_of_zeros_is_no_value(group: dict[str, float]) -> None:
    assert _ssi_labels.label_items(group) == {}


# One item of a real value can be 0
@pytest.mark.parametrize('group', _ZERO_GROUPS, ids=_ZERO_GROUP_IDS)
def test_label_items_keep_a_group_with_one_zero(group: dict[str, float]) -> None:
    label = group | {next(iter(group)): -6.405}
    assert _ssi_labels.label_items(label) == label


# A source that gives two items of a group, both 0, as an index does for the azimuths
def test_label_items_group_of_two_zeros_is_no_value() -> None:
    assert _ssi_labels.label_items({'NORTH_AZIMUTH': 0., 'SUB_SOLAR_AZIMUTH': 0.}) == {}


def test_label_items_keep_a_lone_zero_of_a_group() -> None:
    assert _ssi_labels.label_items({'NORTH_AZIMUTH': 0.}) == {'NORTH_AZIMUTH': 0.}


# DATASET.CAT documents a time from closest approach of "-73T..." as not available
@pytest.mark.parametrize('key', ['TIME_FROM_CLOSEST_APPROACH',
                                 'SATELLITE_TIME_FROM_CLST_APR'])
def test_label_items_unknown_time_from_closest_approach_is_no_value(key: str) -> None:
    assert _ssi_labels.label_items({key: '-73T69:12:71Z'}) == {}


def test_label_items_keep_a_negative_time_from_closest_approach() -> None:
    label = {'NTV_TIME_FROM_CLOSEST_APPROACH': '-073T01:09:57Z'}
    assert _ssi_labels.label_items(label) == label


@pytest.mark.parametrize(('key', 'name'), [
    ('TIME_FROM_CLOSEST_APPROACH', 'NTV_TIME_FROM_CLOSEST_APPROACH'),
    ('SATELLITE_TIME_FROM_CLST_APR', 'NTV_SAT_TIME_FROM_CLOSEST_APRH'),
])
def test_label_items_name_the_time_from_closest_approach(key: str, name: str) -> None:
    assert _ssi_labels.label_items({key: '003T01:09:57Z'}) == {name: '003T01:09:57Z'}

##########################################################################################
# header_items()
##########################################################################################


def test_header_items_strip_strings() -> None:
    items = _ssi_labels.header_items({'ENCODING_TYPE': 'INTEGER COSINE TRANSFORM '})
    assert items['COMPRESSION_TYPE'] == 'INTEGER COSINE TRANSFORM'


@pytest.mark.parametrize('key', ['PARTITION', 'PARTITIO'])
def test_header_items_partition(key: str) -> None:
    items = _ssi_labels.header_items({key: 1})
    assert items == {'SPACECRAFT_CLOCK_CNT_PARTITION': 1}


@pytest.mark.parametrize(('header', 'expected'), [
    ({'EXP': 12.5}, {'EXPOSURE_DURATION': 12.5}),
    ({'PA': 'G1GSGLOBAL02'}, {'OBSERVATION_ID': 'G1GSGLOBAL02', 'ORBIT_NUMBER': 1}),
    ({'TARGET': 'IO'}, {'TARGET_NAME': 'IO'}),
    ({'PICNO': 'G1G0002'}, {'IMAGE_ID': 'G1G0002'}),
    ({'FILTER': 3}, {'FILTER_NUMBER': 3, 'FILTER_NAME': 'VIOLET'}),
    ({'RATE': 5}, {'FRAME_DURATION': 15.166}),
    ({'GAIN': 1}, {'GAIN_MODE_ID': '400K'}),
    ({'PLRANGE': 1.e6}, {'CENTRAL_BODY_DISTANCE': 1.e6}),
], ids=['EXP', 'PA', 'TARGET', 'PICNO', 'FILTER', 'RATE', 'GAIN', 'PLRANGE'])
def test_header_items_give_a_converted_value_only_under_its_pds3_name(
        header: dict[str, Any], expected: dict[str, Any]) -> None:
    assert _ssi_labels.header_items(header) == expected


def test_header_items_clock_count() -> None:
    items = _ssi_labels.header_items({'RIM': 3496321, 'MOD91': 37, 'MOD10': 4, 'MOD8': 6})
    assert items == {'SPACECRAFT_CLOCK_START_COUNT': '03496321.37', 'MOD10': 4,
                     'MOD8': 6}


@pytest.mark.parametrize('header', [{'RATE': 9}, {'GAIN': 0}, {'RIM': 1}],
                         ids=['rate', 'gain', 'mod91'])
def test_header_items_keep_an_unconverted_code(header: dict[str, Any]) -> None:
    assert _ssi_labels.header_items(header) == header


def test_header_items_unknown_filter_number() -> None:
    assert _ssi_labels.header_items({'FILTER': 8}) == {'FILTER_NUMBER': 8}


@pytest.mark.parametrize(('header', 'key'), [
    ({'RATE': 9}, 'FRAME_DURATION'),
    ({'GAIN': 0}, 'GAIN_MODE_ID'),
    ({'FILTER': 8}, 'FILTER_NAME'),
    ({'RIM': 1}, 'SPACECRAFT_CLOCK_START_COUNT'),
    ({'FIBE': ''}, 'EXPOSURE_TYPE'),
    ({'SCETYEAR': 1996}, 'IMAGE_TIME'),
], ids=['rate', 'gain', 'filter', 'mod91', 'flags', 'time-fields'])
def test_header_items_without_a_value(header: dict[str, Any], key: str) -> None:
    assert key not in _ssi_labels.header_items(header)


@pytest.mark.parametrize(('key', 'fill'), [('PICNO', '?'), ('SPICE_C_ID', '?'),
                                           ('TCA', 'CLOSEST TIME'), ('SMEAR', 0.1),
                                           ('RAD', 0.), ('CONE', 0.), ('SMRAZ', 0.),
                                           ('SOLRANGE', 777909100.), ('HSCL', 0.),
                                           ('VSCL', 0.), ('SLRANGE', 0.),
                                           ('PLRANGE', 0.)])
def test_header_items_fill_is_no_value(key: str, fill: Any) -> None:
    assert _ssi_labels.header_items({key: fill}) == {}


@pytest.mark.parametrize('header', [{'RA': 0., 'DEC': 0., 'TWIST': 0.},
                                    {'LAT': 0., 'LON': 0.},
                                    {'NORAZ': 0., 'SUNAZ': 0., 'SCAZ': 0.}],
                         ids=['pointing', 'center', 'azimuths'])
def test_header_items_group_of_zeros_is_no_value(header: dict[str, float]) -> None:
    assert _ssi_labels.header_items(header) == {}


def test_header_items_keep_a_pointing_with_one_zero() -> None:
    items = _ssi_labels.header_items({'RA': -87.6903, 'DEC': 0., 'TWIST': -104.19})
    assert items == {'RIGHT_ASCENSION': -87.6903, 'DECLINATION': 0.,
                     'TWIST_ANGLE': -104.19}


# A fill value of one keyword is a real value of another
def test_header_items_keep_another_keywords_fill() -> None:
    assert _ssi_labels.header_items({'NOTE': 'CLOSEST TIME'}) == {'NOTE': 'CLOSEST TIME'}


def test_header_items_time_from_closest_approach() -> None:
    items = _ssi_labels.header_items({'TCA': '-000T02:21:45Z'})
    assert items == {'NTV_TIME_FROM_CLOSEST_APPROACH': '-000T02:21:45Z'}


def test_header_items_keep_the_flags() -> None:
    assert _ssi_labels.header_items({'FIBE': '1001'}) == {'FIBE': '1001',
                                                          'EXPOSURE_TYPE': 'EXTENDED'}


_SCET = {'SCETYEAR': 1996, 'SCETDAY': 178, 'SCETHOUR': 8, 'SCETMIN': 46, 'SCETSEC': 4,
         'SCETMSEC': 890}


def test_header_items_image_time() -> None:
    assert _ssi_labels.header_items(_SCET) == {'IMAGE_TIME': '1996-06-26T08:46:04.890'}


def test_header_items_earth_received_time() -> None:
    header = {'ERTYEAR': 1996, 'ERTDAY': 32, 'ERTHOUR': 1, 'ERTMIN': 2, 'ERTSEC': 3,
              'ERTMSEC': 4}
    items = _ssi_labels.header_items(header)
    assert items == {'EARTH_RECEIVED_TIME': '1996-02-01T01:02:03.004'}


# A leap second is a valid time
def test_header_items_leap_second() -> None:
    header = _SCET | {'SCETYEAR': 1997, 'SCETDAY': 181, 'SCETHOUR': 23, 'SCETMIN': 59,
                      'SCETSEC': 60}
    assert _ssi_labels.header_items(header)['IMAGE_TIME'] == '1997-06-30T23:59:60.890'


@pytest.mark.parametrize('field', [('SCETYEAR', 0), ('SCETYEAR', -32768),
                                   ('SCETDAY', 367), ('SCETDAY', 0), ('SCETHOUR', 24),
                                   ('SCETMIN', 60), ('SCETSEC', 61), ('SCETMSEC', 1000)],
                         ids=str)
def test_header_items_time_fill_is_no_value(field: tuple[str, int]) -> None:
    assert _ssi_labels.header_items(_SCET | dict([field])) == {}


@pytest.mark.parametrize(('year', 'day'), [(1997, 366), (9999, 366), (10000, 1)],
                         ids=['not-leap', 'last-year', 'beyond-dates'])
def test_header_items_day_beyond_the_year_is_no_value(year: int, day: int) -> None:
    assert _ssi_labels.header_items(_SCET | {'SCETYEAR': year, 'SCETDAY': day}) == {}


def test_header_items_leap_day() -> None:
    header = _SCET | {'SCETYEAR': 1996, 'SCETDAY': 366}
    assert _ssi_labels.header_items(header)['IMAGE_TIME'] == '1996-12-31T08:46:04.890'


def test_header_items_drop_numeric_fill() -> None:
    assert _ssi_labels.header_items({'SMRAZ': -999., 'SMEAR': -999}) == {}

##########################################################################################
# row_items()
##########################################################################################


def test_row_items_give_python_values() -> None:
    row = {'TARGET_NAME': np.str_('IO'), 'FILTER_NUMBER': np.int64(0),
           'ENTROPY': np.float64(2.132), 'CUT_OUT_WINDOW': np.array([1, 1, 400, 400]),
           'TARGET_NAME_mask': np.False_}
    items = _ssi_labels.row_items(row)
    assert json.dumps(items) == ('{"TARGET_NAME": "IO", "FILTER_NUMBER": 0, '
                                 '"ENTROPY": 2.132, "CUT_OUT_WINDOW": [1, 1, 400, 400]}')


# The SSI indexes declare real values as null constants, such as a zero exposure (the
# index's NULL_CONSTANT = 0.), an ON_CHIP_MOSAIC_FLAG of N, and a PRODUCT_VERSION_ID of 0
# (the supplement's NOT_APPLICABLE_CONSTANT); the product labels give the same values
def test_row_items_keep_a_masked_value() -> None:
    row = {'EXPOSURE_DURATION': np.float64(0.), 'EXPOSURE_DURATION_mask': np.True_,
           'ON_CHIP_MOSAIC_FLAG': np.str_('N'), 'ON_CHIP_MOSAIC_FLAG_mask': np.True_,
           'PRODUCT_VERSION_ID': np.int64(0), 'PRODUCT_VERSION_ID_mask': np.True_}
    assert _ssi_labels.row_items(row) == {'EXPOSURE_DURATION': 0.,
                                          'ON_CHIP_MOSAIC_FLAG': 'N',
                                          'PRODUCT_VERSION_ID': 0}


def test_row_items_drop_numeric_fill() -> None:
    assert _ssi_labels.row_items({'SMEAR_AZIMUTH': np.float64(-999.)}) == {}


def test_row_items_drop_nulls() -> None:
    assert _ssi_labels.row_items({'IMAGE_TIME': np.str_('UNK')}) == {}


@pytest.mark.parametrize(('column', 'value'), [
    ('SMEAR_MAGNITUDE', np.float64(0.1)),
    ('SOLAR_DISTANCE', np.float64(777909100.)),
    ('HORIZONTAL_PIXEL_SCALE', np.float64(0.)),
    ('VERTICAL_PIXEL_SCALE', np.float64(0.)),
    ('NTV_TIME_FROM_CLOSEST_APPROACH', np.str_('-73T69:12:71Z ')),
    ('NTV_SAT_TIME_FROM_CLOSEST_APR', np.str_('-73T69:12:71Z ')),
])
def test_row_items_drop_a_placeholder(column: str, value: Any) -> None:
    assert _ssi_labels.row_items({column: value}) == {}


# A corrected table that writes the null constant there gives a real value elsewhere
def test_row_items_keep_a_real_pixel_scale() -> None:
    row = {'HORIZONTAL_PIXEL_SCALE': np.float64(12968.7)}
    assert _ssi_labels.row_items(row) == {'HORIZONTAL_PIXEL_SCALE': 12968.7}


# A value that the table declares as its column's constant is masked by the table reader
@pytest.mark.parametrize(('value', 'mask'), [
    (np.array([-1, -1, -1, -1]), np.array([True, True, True, True])),
    (np.float64(0.), np.True_),
], ids=['window', 'ratio'])
def test_row_items_drop_a_declared_constant(value: Any, mask: Any) -> None:
    assert _ssi_labels.row_items({'X': value, 'X_mask': mask}) == {}


# A corrected table that declares another constant keeps a real value of 0
def test_row_items_keep_an_undeclared_value() -> None:
    row = {'ENCODING_COMPRESSION_RATIO': np.float64(0.),
           'ENCODING_COMPRESSION_RATIO_mask': np.False_}
    assert _ssi_labels.row_items(row) == {'ENCODING_COMPRESSION_RATIO': 0.}


def test_row_items_keep_a_partly_masked_value() -> None:
    row = {'X': np.array([1, -1]), 'X_mask': np.array([False, True])}
    assert _ssi_labels.row_items(row) == {'X': [1, -1]}


@pytest.mark.parametrize('row', [
    {'RIGHT_ASCENSION': np.float64(0.), 'DECLINATION': np.float64(0.),
     'TWIST_ANGLE': np.float64(0.)},
    {'NORTH_AZIMUTH': np.float64(0.), 'SUB_SOLAR_AZIMUTH': np.float64(0.)},
], ids=['pointing', 'azimuths'])
def test_row_items_group_of_zeros_is_no_value(row: dict[str, Any]) -> None:
    assert _ssi_labels.row_items(row) == {}


# Only the declared value that is real is kept; another constant that a corrected table
# declares in the same column is dropped
def test_row_items_drop_another_constant_of_a_real_value_column() -> None:
    row = {'PRODUCT_VERSION_ID': np.int64(-1), 'PRODUCT_VERSION_ID_mask': np.True_}
    assert _ssi_labels.row_items(row) == {}


def test_row_items_name_the_satellite_time_as_the_labels_do() -> None:
    row = {'NTV_SAT_TIME_FROM_CLOSEST_APR': np.str_('000T00:20:02Z')}
    items = _ssi_labels.row_items(row)
    assert items == {'NTV_SAT_TIME_FROM_CLOSEST_APRH': '000T00:20:02Z'}


@pytest.mark.parametrize(('index', 'label'), [
    ('GO-J/JSA-SSI-2-REDR-V1.0', 'GO-J/JSA-SSI-2-REDR-V1.1'),
    ('GO-J/JSA-SSI-4-REDR-V1.0', 'GO-J/JSA-SSI-4-REDR-V1.1'),
    ('GO-A/E-SSI-2-REDR-V1.1', 'GO-A/E-SSI-2-REDR-V1.1'),
], ids=['tour', 'tour-repaired', 'cruise'])
def test_row_items_give_the_labels_data_set_version(index: str, label: str) -> None:
    items = _ssi_labels.row_items({'DATA_SET_ID': np.str_(index.ljust(28))})
    assert items == {'DATA_SET_ID': label}


# The tour indexes write the ratio, minimum, and maximum of GO_0021/C21/IO/C0506406074R
# (3.758, 3.03, 5.428 in its label) as 5.428, 3.758, 3.03
_ROTATED_ROW = {'ENCODING_COMPRESSION_RATIO': np.float64(5.428),
                'ENCODING_MIN_COMPRESSION_RATIO': np.float64(3.758),
                'ENCODING_MAX_COMPRESSION_RATIO': np.float64(3.03)}


def test_row_items_name_the_rotated_compression_ratios() -> None:
    assert _ssi_labels.row_items(_ROTATED_ROW, rotated_ratios=True) == {
        'ENCODING_COMPRESSION_RATIO': 3.758, 'ENCODING_MIN_COMPRESSION_RATIO': 3.03,
        'ENCODING_MAX_COMPRESSION_RATIO': 5.428}


# A corrected table gives each ratio under its own name
def test_row_items_keep_the_ratios_of_a_corrected_table() -> None:
    assert _ssi_labels.row_items(_ROTATED_ROW) == {
        'ENCODING_COMPRESSION_RATIO': 5.428, 'ENCODING_MIN_COMPRESSION_RATIO': 3.758,
        'ENCODING_MAX_COMPRESSION_RATIO': 3.03}


def index_label(order: tuple[str, ...], created: Any = None) -> dict[str, Any]:
    """An index label dictionary whose table has the three ratio columns in the given
    order of their bytes, made at `created` (no PRODUCT_CREATION_TIME if None)."""

    table = {'OBJECT': 'IMAGE_INDEX_TABLE', 'VOLUME_ID': {'START_BYTE': 2}}
    table |= {name: {'NAME': name, 'START_BYTE': 600 + 9 * k}
              for (k, name) in enumerate(order)}
    label: dict[str, Any] = {'^IMAGE_INDEX_TABLE': 'GO_0017_index.tab'}
    if created is not None:
        label['PRODUCT_CREATION_TIME'] = created
    return label | {'IMAGE_INDEX_TABLE': table}


_MIN = 'ENCODING_MIN_COMPRESSION_RATIO'
_MAX = 'ENCODING_MAX_COMPRESSION_RATIO'
_RATIO = 'ENCODING_COMPRESSION_RATIO'


# The column order alone decides: a label made again later with the same order, as every
# regeneration of the metadata tree has been, still has the defect
@pytest.mark.parametrize(('label', 'rotated'), [
    (index_label((_MIN, _MAX, _RATIO)), True),
    (index_label((_MIN, _MAX, _RATIO), dt.datetime(2026, 12, 1, 0, 0, 0)), True),
    (index_label((_RATIO, _MIN, _MAX), dt.datetime(2025, 10, 30, 13, 48, 44)), False),
    (index_label((_MIN, _RATIO, _MAX)), False),
    ({'IMAGE_INDEX_TABLE': {'VOLUME_ID': {'START_BYTE': 2}}}, False),
], ids=['rotated', 'made-later', 'renamed', 'other-order', 'no-ratios'])
def test_ratios_rotated(label: dict[str, Any], rotated: bool) -> None:
    assert _ssi_labels.ratios_rotated(label) is rotated


@pytest.mark.parametrize(('orbit', 'kept'), [(0, True), (35, True), (-99, False)])
def test_row_items_keep_only_tour_orbits(orbit: int, kept: bool) -> None:
    row = {'ORBIT_NUMBER': np.int64(orbit)}
    assert ('ORBIT_NUMBER' in _ssi_labels.row_items(row)) is kept

##########################################################################################
# merge_items() and resolve()
##########################################################################################


def test_merge_label_only() -> None:
    assert _ssi_labels.merge_items({'A': 1}, None) == {'A': 1}


def test_merge_header_only() -> None:
    assert _ssi_labels.merge_items(None, {'A': 1}) == {'A': 1}


def test_merge_keeps_every_key() -> None:
    assert _ssi_labels.merge_items({'A': 1, 'B': 2}, {'B': 2, 'C': 3}) == {'A': 1, 'B': 2,
                                                                           'C': 3}


# A value missing from the winning source is not a disagreement
def test_merge_takes_a_header_item_from_the_label_without_one() -> None:
    merged = _ssi_labels.merge_items({'EXPOSURE_DURATION': 45.83}, {})
    assert merged['EXPOSURE_DURATION'] == 45.83


@pytest.mark.parametrize(('key', 'label', 'header', 'expected'), [
    ('TARGET_NAME', 'J RINGS', 'JUPITER', 'J RINGS'),
    ('CENTRAL_BODY_DISTANCE', 1689682., 1689750., 1689682.),
    ('COMPRESSION_TYPE', 'BARC RATE CONTROL', 'BARC INFORMATION PRESERVING',
     'BARC RATE CONTROL'),
    ('EXPOSURE_DURATION', 45.83, 45.8333, 45.8333),
    ('EXPOSURE_DURATION', 12.5, 12.5003, 12.5003),
    ('TARGET_CENTER_DISTANCE', 665795.8, 665796., 665796.),
    ('X', 3, 4, 3),
    ('X', 3, 3, 3),
], ids=['label-rule', 'label-rule-number', 'label-rule-text', 'header-rule',
        'header-rule-less-digits', 'header-rule-more-digits', 'default', 'agree'])
def test_resolve(key: str, label: Any, header: Any, expected: Any) -> None:
    assert _ssi_labels.resolve(key, label, header) == expected


_HEADER_ITEMS = {
    'CONE_ANGLE', 'DECLINATION', 'EMISSION_ANGLE', 'ENCODING_COMPRESSION_RATIO',
    'ENCODING_MAX_COMPRESSION_RATIO', 'ENCODING_MIN_COMPRESSION_RATIO', 'ENTROPY',
    'EXPOSURE_DURATION', 'INCIDENCE_ANGLE', 'INTERCEPT_POINT_LATITUDE',
    'INTERCEPT_POINT_LINE', 'INTERCEPT_POINT_LINE_SAMPLE', 'INTERCEPT_POINT_LONGITUDE',
    'LOCAL_HOUR_ANGLE', 'NORTH_AZIMUTH', 'PHASE_ANGLE', 'RIGHT_ASCENSION',
    'SUB_SOLAR_AZIMUTH', 'SUB_SOLAR_LATITUDE', 'SUB_SOLAR_LONGITUDE',
    'SUB_SPACECRAFT_AZIMUTH', 'SUB_SPACECRAFT_LATITUDE', 'SUB_SPACECRAFT_LINE',
    'SUB_SPACECRAFT_LINE_SAMPLE', 'SUB_SPACECRAFT_LONGITUDE', 'TARGET_CENTER_DISTANCE',
    'TWIST_ANGLE',
}


def test_the_header_wins_exactly_the_header_items() -> None:
    assert {key for (key, source) in _ssi_labels._RULES.items()
            if source == _ssi_labels._HEADER} == _HEADER_ITEMS

##########################################################################################
