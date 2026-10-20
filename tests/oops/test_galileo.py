##########################################################################################
# tests/oops/test_galileo.py
##########################################################################################
"""Tests for the Galileo mission support, host_galileo.oops._galileo."""

from typing import Any

import julian
import numpy as np
import oops
import pytest
from conftest import (RIM_SECONDS, SCLK_EPOCH, FakeKernel, FakeSpice, FakeSpiceDB,
                      fake_scs2e)
from oops.body import Body

import host_galileo.oops._galileo as galileo_module
from host_galileo.oops._galileo import TIMELINE, _Galileo

##########################################################################################
# Mission timeline
##########################################################################################


def test_timeline_times_follow_their_dates() -> None:
    assert [phase['ET'] for phase in TIMELINE] == pytest.approx(
        [julian.tdb_from_tai(julian.tai_from_iso(phase['UTC'])) for phase in TIMELINE])


def test_timeline_is_in_order() -> None:
    times = [phase['ET'] for phase in TIMELINE]
    assert times == sorted(times)


def test_months_span_the_mission() -> None:
    end = _Galileo.TDB0 + _Galileo.MONTHS * _Galileo.DTDB
    assert end == pytest.approx(_Galileo.TDB1)

##########################################################################################
# initialize() and reset()
##########################################################################################


def test_initialize_defines_the_solar_system_and_spacecraft(fake_spice: FakeSpice) -> None:
    _Galileo.initialize(planets=[5], asof='2000-01-01', irregulars=False)
    assert fake_spice.solar_system == [(('1989-10-18', '2003-09-21'),
                                        {'asof': '2000-01-01', 'planets': [5],
                                         'irregulars': False})]
    assert fake_spice.spice_paths == [('GLL', 'JUPITER')]
    assert fake_spice.leap_seconds == 1


def test_initialize_marks_every_month_loaded() -> None:
    _Galileo.initialize()
    assert _Galileo.CK_LOADED.all()
    assert _Galileo.SPK_LOADED.all()
    assert all(kernels == [] for kernels in _Galileo.CK_LIST)
    assert all(kernels == [] for kernels in _Galileo.SPK_LIST)


def test_initialize_runs_once(fake_spice: FakeSpice) -> None:
    _Galileo.initialize()
    _Galileo.initialize()
    assert len(fake_spice.solar_system) == 1
    assert _Galileo.initialized


def test_reset_clears_state() -> None:
    _Galileo.initialize()
    _Galileo.load_sclk()
    _Galileo.loaded_instruments = ['SSI']
    _Galileo.reset()
    assert not _Galileo.initialized
    assert not _Galileo.sclk_loaded
    assert _Galileo.loaded_instruments == []
    assert not _Galileo.CK_LOADED.any()
    assert not _Galileo.SPK_LOADED.any()

##########################################################################################
# load_kernels(), load_sclk() and tdb_from_sclk()
##########################################################################################


def test_load_kernels_furnishes_the_fixed_set(fake_spice: FakeSpice) -> None:
    _Galileo.load_kernels()
    assert fake_spice.furnished == [_Galileo.SCLK_KERNEL, *galileo_module._KERNELS]


def test_load_sclk_furnishes_the_clock_kernel_once(fake_spice: FakeSpice) -> None:
    _Galileo.load_sclk()
    _Galileo.load_sclk()
    assert fake_spice.furnished == ['Galileo/SCLK/mk00062a.tsc']
    assert fake_spice.leap_seconds == 1


def test_tdb_from_sclk_needs_no_initializer(fake_spice: FakeSpice) -> None:
    assert _Galileo.tdb_from_sclk('00030612.00') == fake_scs2e(-77, '30612:0:0:0')
    assert not _Galileo.initialized
    assert fake_spice.furnished == ['Galileo/SCLK/mk00062a.tsc']


@pytest.mark.parametrize('count', ['00030612:00:0:0', '30612.0', '30612', ' 00030612.00 ',
                                   '30612/0/0/0'],
                         ids=['spice-form', 'short-label', 'rim-only', 'padded',
                              'other-delimiter'])
def test_tdb_from_sclk_equivalent_forms(count: str) -> None:
    assert _Galileo.tdb_from_sclk(count) == _Galileo.tdb_from_sclk('00030612.00')


def test_tdb_from_sclk_passes_every_field() -> None:
    assert _Galileo.tdb_from_sclk('30612.45:3:2') == fake_scs2e(-77, '30612:45:3:2')


def test_tdb_from_sclk_fake_clock_is_linear() -> None:
    assert _Galileo.tdb_from_sclk('1') - _Galileo.tdb_from_sclk('0') == pytest.approx(
        RIM_SECONDS)
    assert _Galileo.tdb_from_sclk('0') == SCLK_EPOCH


@pytest.mark.parametrize('count', ['', 'abc.00', '1:2:3:4:5', '30612.x1'],
                         ids=['empty', 'letters', 'five-fields', 'bad-field'])
def test_tdb_from_sclk_rejects_invalid_counts(count: str) -> None:
    with pytest.raises(oops.OopsValueError,
                       match='invalid Galileo spacecraft clock count') as info:
        _Galileo.tdb_from_sclk(count)
    assert repr(count) in str(info.value)

##########################################################################################
# initialize_kernels()
##########################################################################################


def _month_start(month: float) -> str:
    """The ISO time of a fractional month of the mission, to the millisecond."""

    tdb = _Galileo.TDB0 + month * _Galileo.DTDB
    return str(julian.iso_from_tai(julian.tai_from_tdb(tdb), digits=3))


def test_initialize_kernels_assigns_months() -> None:
    kernel = FakeKernel('a.bc', _month_start(10.5), _month_start(12.5))
    lists = np.empty(_Galileo.MONTHS, dtype='object')
    _Galileo.initialize_kernels([kernel], lists)
    assert [m for m in range(_Galileo.MONTHS) if lists[m]] == [10, 11, 12, 13]


def test_initialize_kernels_clips_to_the_mission() -> None:
    kernel = FakeKernel('a.bc', '1980-01-01T00:00:00', '2010-01-01T00:00:00')
    lists = np.empty(_Galileo.MONTHS, dtype='object')
    _Galileo.initialize_kernels([kernel], lists)
    assert all(kernels == [kernel] for kernels in lists)

##########################################################################################
# load_instruments()
##########################################################################################


def test_load_instruments_includes_ssi_on_the_first_call(fake_spicedb: FakeSpiceDB) -> None:
    _Galileo.load_instruments()
    assert fake_spicedb.named('furnish_inst') == [((-77,), {'inst': ['SSI'], 'asof': None})]
    assert _Galileo.loaded_instruments == ['SSI']


def test_load_instruments_skips_loaded_instruments(fake_spicedb: FakeSpiceDB) -> None:
    _Galileo.load_instruments()
    _Galileo.load_instruments(['SSI'])
    assert len(fake_spicedb.named('furnish_inst')) == 1


def test_load_instruments_leaves_the_callers_list_alone() -> None:
    instruments = ['NIMS']
    _Galileo.load_instruments(instruments)
    assert instruments == ['NIMS']
    assert _Galileo.loaded_instruments == ['NIMS', 'SSI']


def test_load_instruments_formats_asof(fake_spicedb: FakeSpiceDB) -> None:
    _Galileo.load_instruments(asof='2000-02-03')
    assert fake_spicedb.named('furnish_inst')[0][1]['asof'] == '2000-02-03T00:00:00'

##########################################################################################
# spice_instrument_kernel()
##########################################################################################


def test_spice_instrument_kernel(fake_spicedb: FakeSpiceDB) -> None:
    assert _Galileo.spice_instrument_kernel('SSI') == (fake_spicedb.kernel_dict,
                                                       'gll36001.ti')
    assert fake_spicedb.named('select_inst') == [((-77,), {'types': 'IK', 'inst': 'SSI',
                                                           'asof': None})]
    assert fake_spicedb.named('furnish_kernels')[0][1] == {'fast': True}


def test_spice_instrument_kernel_formats_asof(fake_spicedb: FakeSpiceDB) -> None:
    _Galileo.spice_instrument_kernel('SSI', asof='2000-02-03')
    assert fake_spicedb.named('select_inst')[0][1]['asof'] == '2000-02-03T00:00:00'

##########################################################################################
# used_kernels()
##########################################################################################


def _used_bodies(fake_spicedb: FakeSpiceDB, iso: str, **kwargs: Any) -> list[int]:
    """The bodies that used_kernels() asks spicedb about at the given time."""

    tdb = julian.tdb_from_iso(iso)
    assert _Galileo.used_kernels((tdb, tdb + 1.), 'ssi', **kwargs) == fake_spicedb.used
    return list(fake_spicedb.named('used_basenames')[-1][1]['bodies'])


def test_used_kernels_venus_phase(fake_spicedb: FakeSpiceDB) -> None:
    assert _used_bodies(fake_spicedb, '1990-02-10T05:00:00') == [2, 299]


def test_used_kernels_earth_phase_adds_the_moon(fake_spicedb: FakeSpiceDB) -> None:
    assert _used_bodies(fake_spicedb, '1990-12-08T00:00:00') == [3, 399, 301]


def test_used_kernels_asteroid_phase(fake_spicedb: FakeSpiceDB) -> None:
    assert _used_bodies(fake_spicedb, '1993-08-28T00:00:00') == [24310, 2431010]


def test_used_kernels_sl9_phase_has_no_bodies(fake_spicedb: FakeSpiceDB) -> None:
    assert _used_bodies(fake_spicedb, '1994-07-18T00:00:00') == []


def test_used_kernels_jupiter_phase(fake_spicedb: FakeSpiceDB) -> None:
    bodies = _used_bodies(fake_spicedb, '1996-06-26T08:46:04')
    assert bodies == [5, 599, *Body.JUPITER_MOONS_LOADED]


def test_used_kernels_all_planets(fake_spicedb: FakeSpiceDB) -> None:
    bodies = _used_bodies(fake_spicedb, '1996-06-26T08:46:04', return_all_planets=True)
    assert bodies[:6] == [1, 199, 2, 299, 3, 399]
    assert 301 in bodies


def test_used_kernels_passes_the_query(fake_spicedb: FakeSpiceDB) -> None:
    tdb = julian.tdb_from_iso('1996-06-26T08:46:04')
    _Galileo.used_kernels((tdb, tdb + 1.), 'ssi')
    kwargs = fake_spicedb.named('used_basenames')[0][1]
    assert kwargs['types'] is None
    assert kwargs['inst'] == 'ssi'
    assert kwargs['sc'] == -77
    assert kwargs['time'] == (tdb, tdb + 1.)


def test_used_kernels_without_cks(fake_spicedb: FakeSpiceDB) -> None:
    tdb = julian.tdb_from_iso('1996-06-26T08:46:04')
    _Galileo.used_kernels((tdb, tdb + 1.), 'ssi', ck=False)
    types = fake_spicedb.named('used_basenames')[0][1]['types']
    assert 'CK' not in types
    assert 'SPK' in types

##########################################################################################
