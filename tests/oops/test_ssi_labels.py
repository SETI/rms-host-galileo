##########################################################################################
# tests/oops/test_ssi_labels.py
##########################################################################################
"""Tests for the merging of an SSI image's two labels, host_galileo.oops._ssi_labels."""

import datetime as dt
from typing import Any

import pytest

from host_galileo.oops import _ssi_labels


@pytest.mark.parametrize(('value', 'expected'), [
    ('UNK', True), ('N/A', True), ('NULL', True), ('unk', True), (' UNK ', True),
    ('UNKNOWN', False), ('', False), (0, False), (None, False), (-999., False),
])
def test_is_null(value: Any, expected: bool) -> None:
    assert _ssi_labels.is_null(value) is expected


def test_label_items_read_a_set_in_label_order() -> None:
    label = {'CUT_OUT_WINDOW': {2, 1, 399, 392}, 'CUT_OUT_WINDOW_list': [2, 1, 399, 392]}
    assert _ssi_labels.label_items(label)['CUT_OUT_WINDOW'] == [2, 1, 399, 392]


def test_label_items_keep_a_set_without_a_list() -> None:
    assert _ssi_labels.label_items({'X': {1}})['X'] == {1}


def test_label_items_drop_nulls() -> None:
    assert _ssi_labels.label_items({'SMEAR_AZIMUTH': 'UNK', 'X': 1}) == {'X': 1}


@pytest.mark.parametrize(('orbit', 'kept'), [(0, True), (35, True), (36, False),
                                             (12613, False), ('E1', False)])
def test_label_items_keep_only_tour_orbits(orbit: Any, kept: bool) -> None:
    assert ('ORBIT_NUMBER' in _ssi_labels.label_items({'ORBIT_NUMBER': orbit})) is kept


def test_header_items_strip_strings() -> None:
    items = _ssi_labels.header_items({'ENCODING_TYPE': 'INTEGER COSINE TRANSFORM '})
    assert items['COMPRESSION_TYPE'] == 'INTEGER COSINE TRANSFORM'


@pytest.mark.parametrize('key', ['PARTITION', 'PARTITIO'])
def test_header_items_partition(key: str) -> None:
    items = _ssi_labels.header_items({key: 1})
    assert items['SPACECRAFT_CLOCK_CNT_PARTITION'] == 1


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


def test_header_items_earth_received_time() -> None:
    header = {'ERTYEAR': 1996, 'ERTDAY': 32, 'ERTHOUR': 1, 'ERTMIN': 2, 'ERTSEC': 3,
              'ERTMSEC': 4}
    items = _ssi_labels.header_items(header)
    assert items['EARTH_RECEIVED_TIME'] == dt.datetime(1996, 2, 1, 1, 2, 3, 4000)


def test_merge_label_only() -> None:
    assert _ssi_labels.merge_items({'A': 1}, None) == {'A': 1}


def test_merge_header_only() -> None:
    assert _ssi_labels.merge_items(None, {'A': 1}) == {'A': 1}


def test_merge_keeps_every_key() -> None:
    assert _ssi_labels.merge_items({'A': 1, 'B': 2}, {'B': 2, 'C': 3}) == {'A': 1, 'B': 2,
                                                                           'C': 3}


@pytest.mark.parametrize(('label', 'header', 'expected'), [
    (45.83, 45.8333, 45.8333),
    (665795.8, 665796., 665795.8),
    (10.423, 10.4235, 10.4235),
    (0., 5.26301e-06, 5.26301e-06),
    (665800., 665795.8, 665795.8),
    (1689682., 1689750., 1689682.),
    (45.84, 45.8333, 45.84),
    ('J RINGS', 'JUPITER', 'J RINGS'),
    (3, 4, 3),
], ids=['label-rounded', 'header-rounded', 'half-unit', 'zero', 'whole-number',
        'differ', 'not-a-rounding', 'string', 'integer'])
def test_resolve(label: Any, header: Any, expected: Any) -> None:
    assert _ssi_labels.resolve('X', label, header) == expected


def test_resolve_by_a_header_rule(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(_ssi_labels._RULES, 'X', _ssi_labels._HEADER)
    assert _ssi_labels.resolve('X', 'A', 'B') == 'B'


@pytest.mark.parametrize(('number', 'expected'), [
    (45.83, (4, -2)), (665796., (6, 0)), (1689750., (6, 1)), (0., (0, 0)),
    (5.26301e-06, (6, -11)), (-3.25, (3, -2)), (1.5e+20, (2, 19)),
])
def test_precision(number: float, expected: tuple[int, int]) -> None:
    assert _ssi_labels._precision(number) == expected

##########################################################################################
