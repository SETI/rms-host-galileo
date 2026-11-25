##########################################################################################
# tests/oops/test_ssi_from_index.py
##########################################################################################
"""Tests for GalileoSSI.from_index(): the supplement, the pairing of rows, and the
Snapshot of each row."""

import pathlib
from typing import Any

import oops
import pytest
from conftest import (EXPOSURE_MS, INDEX_ROWS, SCLK_COUNT, SSI_PIXELS, SUPPLEMENT_ROWS,
                      fake_scs2e, write_index, write_supplement, write_volume)

from host_galileo.oops.ssi import GalileoSSI

BASENAMES = ['C0346405900R.LBL', 'C0349632100R.LBL', 'C0349632200R.LBL']


def test_from_index_returns_an_iterator(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_volume(tmp_path))
    assert iter(snapshots) is snapshots


def test_from_index_gives_every_row_in_order(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_volume(tmp_path))
    assert [obs.basename for obs in snapshots] == BASENAMES


def test_from_index_keeps_a_zero_exposure(tmp_path: pathlib.Path) -> None:
    snapshots = list(GalileoSSI.from_index(write_volume(tmp_path)))
    assert snapshots[2].texp == 0.


def test_from_index_row_snapshot(tmp_path: pathlib.Path) -> None:
    obs = list(GalileoSSI.from_index(write_volume(tmp_path)))[1]
    assert type(obs) is oops.observation.Snapshot
    assert obs.texp == pytest.approx(EXPOSURE_MS / 1000.)
    assert obs.filter == 'VIOLET'
    assert obs.observation_id == 'G1GSGLOBAL02'
    assert obs.spacecraft_clock_start_count == SCLK_COUNT
    assert obs.window == [1, 0, 400, 392]


def test_from_index_filepath_is_relative(tmp_path: pathlib.Path) -> None:
    obs = list(GalileoSSI.from_index(write_volume(tmp_path)))[1]
    assert obs.filepath == 'GO_0017/G1/GANYMEDE/C0349632100R.LBL'
    assert obs.basename == 'C0349632100R.LBL'


@pytest.mark.parametrize('name', ['image_url', 'image_abspath', 'label_url',
                                  'label_abspath', 'data', 'label_dict', 'header_dict'])
def test_from_index_row_has_no_file_subfield(tmp_path: pathlib.Path, name: str) -> None:
    obs = next(GalileoSSI.from_index(write_volume(tmp_path)))
    assert name not in obs.subfields


def test_from_index_dict_is_the_row_as_the_table_gives_it(tmp_path: pathlib.Path) -> None:
    obs = next(GalileoSSI.from_index(write_volume(tmp_path)))
    assert obs.dict['CUT_OUT_WINDOW'].tolist() == [-1, -1, -1, -1]
    assert obs.dict['CUT_OUT_WINDOW_mask'].tolist() == [True, True, True, True]


def test_from_index_row_without_a_window(tmp_path: pathlib.Path) -> None:
    obs = next(GalileoSSI.from_index(write_volume(tmp_path)))
    assert obs.window == [0, 0, SSI_PIXELS, SSI_PIXELS]
    assert obs.fov is GalileoSSI._FOVS['IM8']


def test_from_index_keeps_the_index_value_of_a_shared_column(
        tmp_path: pathlib.Path) -> None:
    obs = list(GalileoSSI.from_index(write_volume(tmp_path)))[1]
    assert obs.dict['SPACECRAFT_CLOCK_START_COUNT'] == SCLK_COUNT


def test_from_index_row_with_image_time_unk(tmp_path: pathlib.Path) -> None:
    rows = [INDEX_ROWS[0][:7] + ('UNK',) + INDEX_ROWS[0][8:]]
    write_supplement(tmp_path, SUPPLEMENT_ROWS[:1])
    obs = next(GalileoSSI.from_index(write_index(tmp_path, rows)))
    assert obs.time[0] == pytest.approx(fake_scs2e(-77, '3464059:0:0:0') + 1.149,
                                        abs=1.e-6)
    assert obs.time_from_sclk is True


def test_from_index_filter(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_volume(tmp_path),
                                      filter={'FILTER_NAME': 'VIOLET'})
    assert [obs.basename for obs in snapshots] == ['C0349632100R.LBL']


def test_from_index_uses_given_rows(tmp_path: pathlib.Path) -> None:
    rows = [dict(row) for row in oops.Host._index_rows(write_volume(tmp_path),
                                                       validate=())]
    snapshots = GalileoSSI.from_index(tmp_path / 'GO_0017_index.lbl', row_dicts=rows[1:])
    assert [obs.basename for obs in snapshots] == BASENAMES[1:]


def test_from_index_fov_override(tmp_path: pathlib.Path) -> None:
    fov = oops.fov.FlatFOV((1.e-5, 1.e-5), (10, 10))
    obs = next(GalileoSSI.from_index(write_volume(tmp_path), fov=fov))
    assert obs.fov is fov


def test_from_index_frame_override(tmp_path: pathlib.Path) -> None:
    # Unlike the fake scan-platform frame, this one is not an identity matrix
    frame = oops.frame.Cmatrix(oops.Matrix3([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),
                               oops.Frame.J2000, frame_id='GIVEN')
    obs = next(GalileoSSI.from_index(write_volume(tmp_path), frame=frame))
    assert obs.frame is frame


def test_from_index_full_fov(tmp_path: pathlib.Path) -> None:
    obs = list(GalileoSSI.from_index(write_volume(tmp_path), full_fov=True))[1]
    assert obs.fov is GalileoSSI._FOVS['HIS']


def test_from_index_select(tmp_path: pathlib.Path) -> None:
    chosen = list(GalileoSSI.from_index(write_volume(tmp_path), select=slice(None)))
    assert [len(choice) for choice in chosen] == [1, 1, 1]
    assert [choice[0].basename for choice in chosen] == BASENAMES


def test_from_index_given_rows_need_no_supplement(tmp_path: pathlib.Path) -> None:
    rows = [dict(row) for row in oops.Host._index_rows(write_volume(tmp_path),
                                                       validate=())]
    snapshots = GalileoSSI.from_index(tmp_path / 'GO_0017_index.lbl', row_dicts=rows,
                                      supplement=None)
    assert [obs.basename for obs in snapshots] == BASENAMES


def test_from_index_parses_the_labels_with_the_strict_method(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    options: list[Any] = []
    index_rows = oops.Host._index_rows.__func__

    def recording(cls: type, /, filepath: Any, **kwargs: Any) -> Any:
        options.append(kwargs['pds3_options'])
        return index_rows(cls, filepath, **kwargs)

    monkeypatch.setattr(GalileoSSI, '_index_rows', classmethod(recording))
    list(GalileoSSI.from_index(write_volume(tmp_path)))
    assert options == [{'method': 'strict'}]


def test_host_from_index_dispatches_to_ssi(tmp_path: pathlib.Path) -> None:
    snapshots = oops.Host.from_index(write_volume(tmp_path))
    assert [obs.basename for obs in snapshots] == BASENAMES

##########################################################################################
# The supplement
##########################################################################################


def test_from_index_needs_the_supplement_beside_it(tmp_path: pathlib.Path) -> None:
    with pytest.raises(FileNotFoundError):
        GalileoSSI.from_index(write_index(tmp_path))


def test_from_index_without_a_supplement(tmp_path: pathlib.Path) -> None:
    with pytest.raises(oops.host.HostIndexError,
                       match='Galileo SSI index needs its supplemental index'):
        GalileoSSI.from_index(write_volume(tmp_path), supplement=None)


def test_from_index_reads_a_supplement_in_another_directory(
        tmp_path: pathlib.Path) -> None:
    other = tmp_path / 'other'
    other.mkdir()
    supplement = write_supplement(other)
    snapshots = GalileoSSI.from_index(write_index(tmp_path), supplement=supplement)
    assert [obs.basename for obs in snapshots] == BASENAMES


def test_from_index_with_a_missing_supplement(tmp_path: pathlib.Path) -> None:
    with pytest.raises(FileNotFoundError):
        GalileoSSI.from_index(write_volume(tmp_path), supplement=tmp_path / 'nope.lbl')


def test_from_index_cannot_name_the_supplement(tmp_path: pathlib.Path) -> None:
    write_supplement(tmp_path)
    path = write_index(tmp_path, name='GO_0017_table')
    with pytest.raises(oops.host.HostIndexError,
                       match=r'supplemental index required for .*GO_0017_table\.lbl'):
        GalileoSSI.from_index(path)


def test_from_index_of_another_host(tmp_path: pathlib.Path) -> None:
    write_supplement(tmp_path)
    path = write_index(tmp_path, INSTRUMENT_HOST_NAME='CASSINI')
    with pytest.raises(oops.host.WrongHostError, match='not a Galileo SSI index'):
        GalileoSSI.from_index(path)


def test_from_index_with_a_supplement_of_other_files(tmp_path: pathlib.Path) -> None:
    rows = [SUPPLEMENT_ROWS[1], SUPPLEMENT_ROWS[0], SUPPLEMENT_ROWS[2]]
    write_supplement(tmp_path, rows)
    snapshots = GalileoSSI.from_index(write_index(tmp_path))
    with pytest.raises(oops.host.HostIndexError, match='paired by position'):
        next(snapshots)


def test_from_index_with_a_short_supplement(tmp_path: pathlib.Path) -> None:
    write_supplement(tmp_path, SUPPLEMENT_ROWS[:2])
    snapshots = GalileoSSI.from_index(write_index(tmp_path))
    with pytest.raises(oops.host.HostIndexError,
                       match='the supplement ends before index row 2'):
        list(snapshots)


def test_from_index_unknown_keyword(tmp_path: pathlib.Path) -> None:
    with pytest.raises(TypeError, match="unexpected keyword argument 'astrometry'"):
        GalileoSSI.from_index(write_volume(tmp_path),
                              astrometry=True)

##########################################################################################
