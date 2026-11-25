##########################################################################################
# tests/oops/test_ssi_labels.py
##########################################################################################
"""Tests for the SSI metadata items, host_galileo.oops._ssi_labels."""

import datetime as dt
from typing import Any

import numpy as np
import pytest
import vicar

from host_galileo.oops import _ssi_labels


def header(*params: tuple[str, Any]) -> vicar.VicarLabel:
    """A VICAR header holding the system keywords of a tiny image, then `params`."""

    label = vicar.VicarImage.from_array(np.zeros((2, 2), dtype='uint8')).label
    for param in params:
        label.append(param)
    return label


_SCET = [('SCETYEAR', 1996), ('SCETDAY', 178), ('SCETHOUR', 8), ('SCETMIN', 46),
         ('SCETSEC', 4), ('SCETMSEC', 890)]

##########################################################################################
# label_items()
##########################################################################################


def test_label_items_keep_every_value() -> None:
    items = _ssi_labels.label_items({'TARGET_NAME': 'IO', 'FILTER_NUMBER': 3,
                                     'IMAGE_TIME_fmt': '1996-06-26T08:46:04.890'})
    assert items == {'TARGET_NAME': 'IO', 'FILTER_NUMBER': 3,
                     'IMAGE_TIME_fmt': '1996-06-26T08:46:04.890'}


@pytest.mark.parametrize('value', [None, 'UNK', 'N/A', 'not applicable', '---', ''],
                         ids=['none', 'unk', 'na', 'not-applicable', 'dashes', 'empty'])
def test_label_items_leave_out_a_null(value: Any) -> None:
    assert 'SMEAR_AZIMUTH' not in _ssi_labels.label_items({'SMEAR_AZIMUTH': value})


@pytest.mark.parametrize('key', ['CUT_OUT_WINDOW', 'TRUTH_WINDOW'])
def test_label_items_give_an_ordered_set_in_label_order(key: str) -> None:
    items = _ssi_labels.label_items({key: {392, 1, 2, 399},
                                     key + '_list': [2, 1, 399, 392]})
    assert items[key] == [2, 1, 399, 392]


@pytest.mark.parametrize('orbit', [0, 1, 35], ids=['J0', 'G1', 'orbit-35'])
def test_label_items_keep_an_orbit_of_the_tour(orbit: int) -> None:
    assert _ssi_labels.label_items({'ORBIT_NUMBER': orbit})['ORBIT_NUMBER'] == orbit


@pytest.mark.parametrize('orbit', [36, 12613, -1], ids=['after-tour', 'E1', 'negative'])
def test_label_items_leave_out_an_orbit_outside_the_tour(orbit: int) -> None:
    assert 'ORBIT_NUMBER' not in _ssi_labels.label_items({'ORBIT_NUMBER': orbit})

##########################################################################################
# header_items()
##########################################################################################


def test_header_items_keep_an_unconverted_parameter() -> None:
    items = _ssi_labels.header_items(header(('SPICE_C_ID', 'M905')))
    assert items['SPICE_C_ID'] == 'M905'


@pytest.mark.parametrize('value', [-999, -999., 'NOT APPLICABLE'],
                         ids=['int-fill', 'real-fill', 'null-text'])
def test_header_items_leave_out_a_fill(value: Any) -> None:
    assert 'SMRAZ' not in _ssi_labels.header_items(header(('SMRAZ', value)))


def test_header_items_give_a_repeated_parameter_once() -> None:
    items = _ssi_labels.header_items(header(('TASK', 'SSIMERGE'), ('TASK', 'CATLABEL')))
    assert items['TASK'] == 'SSIMERGE'


@pytest.mark.parametrize(('vicar_key', 'pds3_key'), [
    ('EXP', 'EXPOSURE_DURATION'),
    ('PA', 'OBSERVATION_ID'),
    ('PARTITION', 'SPACECRAFT_CLOCK_CNT_PARTITION'),
    ('PARTITIO', 'SPACECRAFT_CLOCK_CNT_PARTITION'),
], ids=['EXP', 'PA', 'PARTITION', 'PARTITIO'])
def test_header_items_give_a_pair_under_its_pds3_name(vicar_key: str,
                                                      pds3_key: str) -> None:
    assert _ssi_labels.header_items(header((vicar_key, 7)))[pds3_key] == 7


def test_header_items_drop_the_vicar_name_of_a_pair() -> None:
    assert 'EXP' not in _ssi_labels.header_items(header(('EXP', 262.5)))


def test_header_items_give_the_clock_count() -> None:
    items = _ssi_labels.header_items(header(('RIM', 3496321), ('MOD91', 5)))
    assert items['SPACECRAFT_CLOCK_START_COUNT'] == '03496321.05'


@pytest.mark.parametrize('key', ['RIM', 'MOD91'])
def test_header_items_drop_the_clock_fields(key: str) -> None:
    assert key not in _ssi_labels.header_items(header(('RIM', 3496321), ('MOD91', 5)))


@pytest.mark.parametrize('key', ['MOD10', 'MOD8'])
def test_header_items_keep_the_finer_clock_fields(key: str) -> None:
    assert _ssi_labels.header_items(header((key, 4)))[key] == 4


def test_header_items_name_the_filter() -> None:
    assert _ssi_labels.header_items(header(('FILTER', 3)))['FILTER_NAME'] == 'VIOLET'


def test_header_items_leave_an_unknown_filter_unnamed() -> None:
    assert 'FILTER_NAME' not in _ssi_labels.header_items(header(('FILTER', 9)))


@pytest.mark.parametrize(('vicar_key', 'code', 'pds3_key', 'value'), [
    ('RATE', 5, 'FRAME_DURATION', 15.166),
    ('GAIN', 2, 'GAIN_MODE_ID', '100K'),
], ids=['RATE', 'GAIN'])
def test_header_items_convert_a_code(vicar_key: str, code: int, pds3_key: str,
                                     value: Any) -> None:
    assert _ssi_labels.header_items(header((vicar_key, code)))[pds3_key] == value


def test_header_items_drop_a_converted_code() -> None:
    assert 'RATE' not in _ssi_labels.header_items(header(('RATE', 5)))


def test_header_items_keep_an_unknown_code() -> None:
    assert _ssi_labels.header_items(header(('RATE', 9)))['RATE'] == 9


@pytest.mark.parametrize(('params', 'expected'), [
    ((('FIBE', '0001'),), 'EXTENDED'),
    ((('FIBE', '0000'),), 'NORMAL'),
    ((('FIBE', '0001'), ('MOFIBE', '001000')), 'NORMAL'),
], ids=['FIBE-extended', 'FIBE-normal', 'MOFIBE-first'])
def test_header_items_give_the_exposure_type(params: tuple[tuple[str, str], ...],
                                             expected: str) -> None:
    assert _ssi_labels.header_items(header(*params))['EXPOSURE_TYPE'] == expected


def test_header_items_give_the_image_time() -> None:
    items = _ssi_labels.header_items(header(*_SCET))
    assert items['IMAGE_TIME'] == dt.datetime(1996, 6, 26, 8, 46, 4, 890000)


def test_header_items_give_the_earth_received_time() -> None:
    items = _ssi_labels.header_items(header(('ERTYEAR', 1996), ('ERTDAY', 193),
                                            ('ERTHOUR', 10), ('ERTMIN', 44),
                                            ('ERTSEC', 54), ('ERTMSEC', 307)))
    assert items['EARTH_RECEIVED_TIME'] == dt.datetime(1996, 7, 11, 10, 44, 54, 307000)


def test_header_items_drop_the_time_fields() -> None:
    assert 'SCETDAY' not in _ssi_labels.header_items(header(*_SCET))


def test_header_items_give_no_time_for_the_year_fill() -> None:
    params = [('SCETYEAR', -32768), *_SCET[1:]]
    assert 'IMAGE_TIME' not in _ssi_labels.header_items(header(*params))


def test_header_items_give_no_time_without_every_field() -> None:
    assert 'IMAGE_TIME' not in _ssi_labels.header_items(header(*_SCET[:-1]))


def test_header_items_derive_no_orbit_number() -> None:
    assert 'ORBIT_NUMBER' not in _ssi_labels.header_items(header(('PA', 'G1GSGLOBAL02')))

##########################################################################################
# merge_items()
##########################################################################################


@pytest.mark.parametrize(('key', 'expected'), [
    ('EXPOSURE_DURATION', 'header'),
    ('TARGET_NAME', 'label'),
    ('FILTER_NAME', 'label'),
], ids=['header-rule', 'label-rule', 'no-rule'])
def test_merge_items_keep_the_value_of_the_rule(key: str, expected: str) -> None:
    merged = _ssi_labels.merge_items({key: 'label'}, {key: 'header'})
    assert merged[key] == expected


def test_merge_items_fall_to_the_label_when_the_header_lacks_an_item() -> None:
    merged = _ssi_labels.merge_items({'EXPOSURE_DURATION': 262.5}, {'MOD10': 0})
    assert merged['EXPOSURE_DURATION'] == 262.5


def test_merge_items_fall_to_the_header_when_the_label_lacks_an_item() -> None:
    merged = _ssi_labels.merge_items({'FILTER_NAME': 'VIOLET'}, {'TARGET_NAME': 'IO'})
    assert merged['TARGET_NAME'] == 'IO'


@pytest.mark.parametrize(('label', 'header_items', 'expected'), [
    ({'A': 1}, None, {'A': 1}),
    (None, {'B': 2}, {'B': 2}),
], ids=['label-only', 'header-only'])
def test_merge_items_of_one_source(label: dict[str, int] | None,
                                   header_items: dict[str, int] | None,
                                   expected: dict[str, int]) -> None:
    assert _ssi_labels.merge_items(label, header_items) == expected

##########################################################################################
