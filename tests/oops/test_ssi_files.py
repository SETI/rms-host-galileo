##########################################################################################
# tests/oops/test_ssi_files.py
##########################################################################################
"""Tests that the SSI image, its label, and its index row give consistent Snapshots.

The files in ``test_files`` are cut down from GO_0022/I24/IO/C0520792800R: the image ends
after its VICAR header, and each index keeps its first 30 rows. The image therefore gives
no data, so each file is read for astrometry: the label alone, or the VICAR header alone.
Where the two labels agree, the Snapshots agree. The exposure, which the PDS3 label gives
to two decimal places and the VICAR header to four, is the one value that the labels give
differently and that changes a Snapshot attribute.
"""

import pathlib
from typing import Any

import julian
import oops
import pytest
from pdsparser import Pds3Label
from vicar import VicarLabel

from host_galileo.oops import _ssi_labels
from host_galileo.oops.ssi import GalileoSSI

TEST_FILES = pathlib.Path(__file__).parents[2] / 'test_files'
INDEX = TEST_FILES / 'GO_0022_index_edited.lbl'
INDEX_ROWS = 30

SOURCES = ['IMG', 'index']

# The exposure in seconds, as the PDS3 label and the VICAR header give it
LABEL_EXPOSURE = 0.04583
HEADER_EXPOSURE = 0.0458333

# The middle of the exposure in seconds TDB: IMAGE_TIME in the label and the index, and
# the SCET fields in the VICAR header
IMAGE_MIDTIME = float(julian.tdb_from_tai(julian.tai_from_iso('1999-10-11T04:29:52.510')))

# The corners and center of the 800x400 field of view
CORNERS = oops.Pair([(0, 0), (800, 0), (0, 400), (800, 400), (400, 200)])


def load(source: str) -> Any:
    """The Snapshot of C0520792800R from the image ('IMG'), label ('LBL') or 'index'."""

    if source == 'index':
        return next(obs for obs in GalileoSSI.from_index(INDEX)
                    if obs.basename == 'C0520792800R.IMG')
    return GalileoSSI.from_file(TEST_FILES / f'C0520792800R_edited.{source}',
                                astrometry=True)


def items() -> tuple[dict[str, Any], dict[str, Any]]:
    """The metadata items of the test label and of the test image's VICAR header."""

    label = Pds3Label(TEST_FILES / 'C0520792800R_edited.LBL', sources=True).as_dict()
    header = VicarLabel(TEST_FILES / 'C0520792800R_edited.IMG').as_dict(unique=False)
    return (_ssi_labels.label_items(label), _ssi_labels.header_items(header))


def test_index_reads_every_row() -> None:
    assert len(GalileoSSI.from_index(INDEX)) == INDEX_ROWS


def test_index_merges_the_supplement() -> None:
    assert load('index').dict['TELEMETRY_FORMAT_ID'] == 'IM8'


@pytest.mark.parametrize('source', SOURCES)
@pytest.mark.parametrize('attribute', ['instrument', 'filter', 'observation_id',
                                       'label_target', 'target', 'lightsource', 'window',
                                       'shape', 'swap_uv', 'full_fov', 'time_from_sclk',
                                       'spacecraft_clock_start_count',
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


@pytest.mark.parametrize(('source', 'texp'), [('LBL', LABEL_EXPOSURE),
                                              ('index', LABEL_EXPOSURE),
                                              ('IMG', HEADER_EXPOSURE)])
def test_exposure_comes_from_the_label_read(source: str, texp: float) -> None:
    assert load(source).texp == pytest.approx(texp, rel=1.e-12)


@pytest.mark.parametrize(('source', 'texp'), [('LBL', LABEL_EXPOSURE),
                                              ('index', LABEL_EXPOSURE),
                                              ('IMG', HEADER_EXPOSURE)])
def test_exposure_is_centered_on_the_image_time(source: str, texp: float) -> None:
    assert load(source).time[0] == pytest.approx(IMAGE_MIDTIME - texp / 2., abs=1.e-6)


@pytest.mark.parametrize('attribute', ['basename', 'image_url'])
def test_label_names_the_image(attribute: str) -> None:
    assert getattr(load('LBL'), attribute) == getattr(load('IMG'), attribute)


@pytest.mark.parametrize('source', SOURCES)
def test_fov_shape_matches_the_label(source: str) -> None:
    expected = load('LBL').fov.uv_shape.vals.tolist()
    assert load(source).fov.uv_shape.vals.tolist() == expected


@pytest.mark.parametrize('source', SOURCES)
def test_fov_los_matches_the_label(source: str) -> None:
    los = load(source).fov.los_from_uv(CORNERS).vals
    assert los == pytest.approx(load('LBL').fov.los_from_uv(CORNERS).vals, abs=1.e-15)


def test_labels_differ_only_by_rounding() -> None:
    (label, header) = items()
    differ = [key for key in label.keys() & header.keys()
              if label[key] != header[key]
              and not _ssi_labels._is_rounded(label[key], header[key])
              and not _ssi_labels._is_rounded(header[key], label[key])]
    assert differ == []


def test_merge_keeps_the_precise_exposure() -> None:
    assert _ssi_labels.merge_items(*items())['EXPOSURE_DURATION'] == 45.8333


def test_merge_keeps_the_precise_distance() -> None:
    assert _ssi_labels.merge_items(*items())['TARGET_CENTER_DISTANCE'] == 2863.583

##########################################################################################
