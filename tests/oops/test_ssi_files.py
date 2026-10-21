##########################################################################################
# tests/oops/test_ssi_files.py
##########################################################################################
"""Tests that the SSI image, its label, and its index row all give the same Snapshot.

The files in ``test_files`` are cut down from GO_0022/I24/IO/C0520792800R: the image ends
after its VICAR header, and each index keeps its first 30 rows.
"""

import pathlib
from typing import Any

import oops
import pytest

from host_galileo.oops.ssi import GalileoSSI

TEST_FILES = pathlib.Path(__file__).parents[2] / 'test_files'
INDEX = TEST_FILES / 'GO_0022_index_edited.lbl'
INDEX_ROWS = 30

# The VICAR header gives the exposure as 45.8333 ms, the PDS3 label and index as 45.83 ms
EXPOSURE_TOLERANCE = 1.e-5      # seconds

SOURCES = ['IMG', 'index']

# The corners and center of the 800x400 field of view
CORNERS = oops.Pair([(0, 0), (800, 0), (0, 400), (800, 400), (400, 200)])


def load(source: str) -> Any:
    """The Snapshot of C0520792800R from the image ('IMG'), label ('LBL') or 'index'."""

    if source == 'index':
        return next(obs for obs in GalileoSSI.from_index(INDEX)
                    if obs.basename == 'C0520792800R.LBL')
    return GalileoSSI.from_file(TEST_FILES / f'C0520792800R_edited.{source}',
                                astrometry=True)


def test_index_reads_every_row() -> None:
    assert len(GalileoSSI.from_index(INDEX)) == INDEX_ROWS


def test_index_merges_the_supplement() -> None:
    assert load('index').dict['TELEMETRY_FORMAT_ID'] == 'IM8'


@pytest.mark.parametrize('source', SOURCES)
@pytest.mark.parametrize('attribute', ['instrument', 'filter', 'observation_id',
                                       'label_target', 'target', 'lightsource', 'window',
                                       'shape', 'swap_uv', 'full_fov', 'time_from_sclk',
                                       'spice_frame_name', 'spice_kernels',
                                       'calibrations', 'u_axis', 'v_axis', 't_axis'])
def test_attribute_matches_the_label(source: str, attribute: str) -> None:
    assert getattr(load(source), attribute) == getattr(load('LBL'), attribute)


@pytest.mark.parametrize('source', SOURCES)
def test_path_matches_the_label(source: str) -> None:
    assert load(source).path.path_id == load('LBL').path.path_id


@pytest.mark.parametrize('source', SOURCES)
def test_frame_matches_the_label(source: str) -> None:
    assert load(source).frame.frame_id == load('LBL').frame.frame_id


@pytest.mark.parametrize('source', SOURCES)
def test_midtime_matches_the_label(source: str) -> None:
    assert load(source).midtime == pytest.approx(load('LBL').midtime, abs=1.e-6)


@pytest.mark.parametrize('source', SOURCES)
def test_start_time_matches_the_label(source: str) -> None:
    assert load(source).time[0] == pytest.approx(load('LBL').time[0],
                                                 abs=EXPOSURE_TOLERANCE)


@pytest.mark.parametrize('source', SOURCES)
def test_exposure_matches_the_label(source: str) -> None:
    assert load(source).texp == pytest.approx(load('LBL').texp, abs=EXPOSURE_TOLERANCE)


@pytest.mark.parametrize('source', SOURCES)
def test_fov_shape_matches_the_label(source: str) -> None:
    assert load(source).fov.uv_shape.vals.tolist() == load('LBL').fov.uv_shape.vals.tolist()


@pytest.mark.parametrize('source', SOURCES)
def test_fov_los_matches_the_label(source: str) -> None:
    los = load(source).fov.los_from_uv(CORNERS).vals
    assert los == pytest.approx(load('LBL').fov.los_from_uv(CORNERS).vals, abs=1.e-15)

##########################################################################################
