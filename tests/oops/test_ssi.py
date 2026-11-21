##########################################################################################
# tests/oops/test_ssi.py
##########################################################################################
"""Tests for the Galileo SSI host, host_galileo.oops.ssi."""

from collections.abc import Iterator
import itertools
import json
import os
import pathlib
import re
import subprocess
import sys
from typing import Any, cast

import cspyce
import julian
import numpy as np
import oops
from filecache import FCPath
from oops.host._index_rows import CHUNK_ROWS
from pdsparser import Pds3Label
import pytest
import vicar
from conftest import (SSI_DISTORTION, SSI_FOCAL_LENGTH, SSI_PIXEL_SIZE, SSI_PIXELS,
                      FakeSpice, FakeSpiceDB, fake_scs2e, pds3_lines, write_text)

import host_galileo
import host_galileo.oops
import host_galileo.oops.ssi as ssi
from host_galileo.oops          import _ssi_labels
from host_galileo.oops._galileo import _Galileo
from host_galileo.oops.ssi      import GalileoSSI

IMAGE_TIME = '1996-06-26T08:46:04.890Z'
SCLK_COUNT = '03496321.00'
EXPOSURE_MS = 262.5

# The label of GO_0017/G1/GANYMEDE/C0349632100R, reduced to the keywords the host reads
_LABEL = {
    'SPACECRAFT_NAME'             : 'GALILEO ORBITER',
    'INSTRUMENT_NAME'             : 'SOLID STATE IMAGING SYSTEM',
    'SPACECRAFT_CLOCK_START_COUNT': SCLK_COUNT,
    'IMAGE_TIME'                  : IMAGE_TIME,
    'OBSERVATION_ID'              : 'G1GSGLOBAL02',
    'TARGET_NAME'                 : 'GANYMEDE',
    'FILTER_NAME'                 : 'VIOLET',
    'EXPOSURE_DURATION'           : EXPOSURE_MS,
    'FRAME_DURATION'              : 15.166,
    'EXPOSURE_TYPE'               : 'NORMAL',
    'TELEMETRY_FORMAT_ID'         : 'HIS',
    'CUT_OUT_WINDOW'              : [2, 1, 399, 392],
}


def label_dict(**changes: Any) -> dict[str, Any]:
    """The test label as a dictionary; a value of None removes the keyword."""

    dict_ = {**_LABEL, **changes}
    return {key: value for (key, value) in dict_.items() if value is not None}


def image_array() -> np.ndarray:
    """The 800x800 test image, whose values identify their pixels."""

    return (np.arange(SSI_PIXELS**2, dtype='uint32').reshape(SSI_PIXELS, SSI_PIXELS)
            % 251).astype('uint8')


def write_pds3_label(directory: pathlib.Path, stem: str, *, record: int = 1,
                     samples: int = SSI_PIXELS, extra_lines: tuple[str, ...] = (),
                     **changes: Any) -> pathlib.Path:
    """Write a detached PDS3 label for the 800x800 image `stem`.IMG; return its path.

    The label says that the image begins at `record` of the file, counted from one, and
    that each line has `samples` samples. Any `extra_lines` are added as written.
    """

    lines = ['PDS_VERSION_ID = PDS3', 'RECORD_TYPE = FIXED_LENGTH',
             f'RECORD_BYTES = {samples}', f'FILE_RECORDS = {SSI_PIXELS}',
             f'^IMAGE = ("{stem}.IMG", {record})', *extra_lines]
    lines += pds3_lines(label_dict(**changes))
    lines += ['OBJECT = IMAGE', f'  LINES = {SSI_PIXELS}',
              f'  LINE_SAMPLES = {samples}', '  SAMPLE_BITS = 8',
              '  SAMPLE_TYPE = UNSIGNED_INTEGER', 'END_OBJECT = IMAGE']
    return write_text(directory / f'{stem}.LBL', lines)


# HOST is the archive's, not the machine's that writes or reads the test file
_VICAR_LABEL: dict[str, int | float | str | list[int]] = {
    'HOST'          : 'AXP-VMS',
    'MISSION'       : 'GALILEO',
    'SENSOR'        : 'SSI',
    'RIM'           : 3496321,
    'MOD91'         : 0,
    'MOD10'         : 0,
    'MOD8'          : 0,
    'PA'            : 'G1GSGLOBAL02',
    'TARGET'        : 'GANYMEDE',
    'FILTER'        : 3,
    'EXP'           : EXPOSURE_MS,
    'RATE'          : 5,
    'GAIN'          : 1,
    'PICNO'         : 'G1G0002',
    'MOFIBE'        : '001000',
    'TLMFMT'        : 'HIS',
    'CUT_OUT_WINDOW': [2, 1, 399, 392],
    'SCETYEAR'      : 1996,
    'SCETDAY'       : 178,
    'SCETHOUR'      : 8,
    'SCETMIN'       : 46,
    'SCETSEC'       : 4,
    'SCETMSEC'      : 890,
}


def write_vicar_image(directory: pathlib.Path, name: str = 'C0349632100R.IMG', /,
                      **extra: float | str | None) -> pathlib.Path:
    """Write an image as a VICAR file `name` with no PDS3 label; return the image path.

    The header holds the VICAR keywords of C0349632100R that the host reads, changed or
    extended by any `extra` keywords; a value of None removes the keyword.
    """

    vic = vicar.VicarImage.from_array(image_array())
    for (key, value) in {**_VICAR_LABEL, **extra}.items():
        if value is not None:
            vic[key] = value
    path = directory / name
    vic.write_file(path)
    return path


def write_pds3_image(directory: pathlib.Path, *,
                     header: dict[str, float | str | None] | None = None,
                     image_name: str = 'C0349632100R.IMG', record_error: int = 0,
                     samples: int = SSI_PIXELS, extra_lines: tuple[str, ...] = (),
                     **changes: Any) -> pathlib.Path:
    """Write an image as a VICAR file with its detached PDS3 label; return the label path.

    The VICAR header is changed by `header` as :func:`write_vicar_image` changes it by
    its keywords, and the image file is named `image_name`, whatever case the label uses.
    The label's IMAGE pointer names the first record after the VICAR header, as an archive
    label's does, plus `record_error`. The label's other keywords are changed by
    `changes` as :func:`label_dict` changes them.
    """

    image = write_vicar_image(directory, image_name, **(header or {}))
    header_bytes = cast(int, vicar.VicarLabel(image)['LBLSIZE'])
    record = header_bytes // SSI_PIXELS + 1 + record_error
    return write_pds3_label(directory, 'C0349632100R', record=record, samples=samples,
                            extra_lines=extra_lines, **changes)


OVERRIDE_NAMES = ('target', 'lightsource', 'path', 'frame', 'cadence', 'fov',
                  'calibrations', 'timeshift', 'navigation', 'tracker', 'parallel')


def make_snapshot(dict_: dict[str, Any], *, filepath: str = 'X.IMG', data: Any = None,
                  full_fov: bool = False, return_all_planets: bool = False,
                  default_mode: str | None = None,
                  **overrides: Any) -> oops.observation.Snapshot:
    """A Snapshot from `_make_snapshot` for the image file `filepath`, with every override
    not given set to None."""

    GalileoSSI.initialize()
    obs: oops.observation.Snapshot = GalileoSSI._make_snapshot(
        dict_, filepath=FCPath(filepath), data=data, full_fov=full_fov,
        return_all_planets=return_all_planets, default_mode=default_mode,
        overrides={name: overrides.get(name) for name in OVERRIDE_NAMES})
    return obs


def snapshot(**changes: Any) -> oops.observation.Snapshot:
    """A Snapshot made from the test label as an index row, without any file."""

    return make_snapshot(_ssi_labels.row_items(label_dict(**changes)),
                         filepath='GO_0017/X.IMG')


def tdb_from_image_time(image_time: str = IMAGE_TIME) -> float:
    """The time in seconds TDB of an IMAGE_TIME."""

    return float(julian.tdb_from_tai(julian.tai_from_iso(image_time)))

##########################################################################################
# Registration and package layout
##########################################################################################


def test_ssi_is_registered_with_oops() -> None:
    assert oops.Host._LOOKUP['Galileo SSI'] is GalileoSSI


def test_module_exports() -> None:
    assert ssi.__all__ == ['GalileoSSI']


@pytest.mark.parametrize('module', [host_galileo, host_galileo.oops],
                         ids=['host_galileo', 'oops'])
def test_package_exposes_host(module: Any) -> None:
    assert module.GalileoSSI is GalileoSSI


# Prints whether oops and host_galileo.oops.ssi are loaded after importing host_galileo,
# and again after its GalileoSSI is used.
_LAZY_IMPORT_SCRIPT = """
import sys
names = ('oops', 'host_galileo.oops.ssi')
import host_galileo
print(*(name in sys.modules for name in names))
from host_galileo import GalileoSSI
print(*(name in sys.modules for name in names))
print(sys.modules['oops'].Host._LOOKUP['Galileo SSI'] is GalileoSSI)
"""


def test_host_is_imported_lazily() -> None:
    src = pathlib.Path(__file__).parents[2] / 'src'
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join([str(src), env.get('PYTHONPATH', '')])
    result = subprocess.run([sys.executable, '-c', _LAZY_IMPORT_SCRIPT], env=env,
                            capture_output=True, text=True, check=True)
    assert result.stdout.splitlines() == ['False False', 'True True', 'True']


def test_unknown_attribute() -> None:
    with pytest.raises(AttributeError,
                       match="module 'host_galileo' has no attribute 'Nope'"):
        _ = host_galileo.Nope

##########################################################################################
# Detectors
##########################################################################################


@pytest.mark.parametrize(('label', 'expected'), [
    ({'SPACECRAFT_NAME': 'GALILEO ORBITER',
      'INSTRUMENT_NAME': 'SOLID STATE IMAGING SYSTEM'}, True),
    ({'SPACECRAFT_NAME': 'GALILEO ORBITER',
      'INSTRUMENT_NAME': 'SOLID_STATE_IMAGING'},        True),
    ({'SPACECRAFT_NAME': 'GALILEO ORBITER',
      'INSTRUMENT_NAME': 'NEAR INFRARED MAPPING'},      False),
    ({'SPACECRAFT_NAME': 'CASSINI ORBITER',
      'INSTRUMENT_NAME': 'SOLID STATE IMAGING'},        False),
    ({}, False),
], ids=['spaces', 'underscores', 'other-instrument', 'other-spacecraft', 'empty'])
def test_detect_in_pds3(label: dict[str, str], expected: bool) -> None:
    assert GalileoSSI._detect_in_pds3(label) is expected


@pytest.mark.parametrize(('label', 'expected'), [
    ({'MISSION': 'GALILEO', 'SENSOR': 'SSI'}, True),
    ({'MISSION': 'GALILEO', 'SENSOR': 'NIMS'}, False),
    ({'MISSION': 'CASSINI', 'SENSOR': 'SSI'}, False),
], ids=['ssi', 'other-sensor', 'other-mission'])
def test_detect_in_vicar(label: dict[str, str], expected: bool) -> None:
    assert GalileoSSI._detect_in_vicar(label) is expected


@pytest.mark.parametrize(('label', 'expected'), [
    ({'INSTRUMENT_HOST_NAME': 'GALILEO', 'INSTRUMENT_ID': 'SSI',
      'IMAGE_INDEX_TABLE': {'ROWS': 1}}, True),
    ({'INSTRUMENT_HOST_NAME': 'GALILEO', 'INSTRUMENT_ID': 'NIMS'}, False),
    ({'INSTRUMENT_HOST_NAME': 'CASSINI ORBITER', 'INSTRUMENT_ID': 'SSI'}, False),
    ({'IMAGE_INDEX_TABLE': {'ROWS': 1}}, None),
], ids=['ssi', 'other-instrument', 'other-host', 'unnamed'])
def test_detect_in_index(label: dict[str, Any], expected: bool | None) -> None:
    assert GalileoSSI._detect_in_index(label) is expected


@pytest.mark.parametrize(('row', 'expected'), [
    ({'MISSION_NAME': 'GALILEO', 'INSTRUMENT_ID': 'SSI'}, True),
    ({'MISSION_NAME': 'GALILEO', 'INSTRUMENT_ID': 'NIMS'}, False),
    ({}, False),
], ids=['ssi', 'other-instrument', 'empty'])
def test_detect_in_row(row: dict[str, str], expected: bool) -> None:
    assert GalileoSSI._detect_in_row(row) is expected

##########################################################################################
# initialize(), FOVs and frames
##########################################################################################


def test_initialize_loads_the_mission(fake_spice: FakeSpice,
                                      fake_spicedb: FakeSpiceDB) -> None:
    GalileoSSI.initialize(planets=[5], asof='2000-01-01')
    assert fake_spice.solar_system[0][1]['planets'] == [5]
    assert fake_spicedb.named('furnish_inst')[0][1]['inst'] == ['SSI']
    assert fake_spicedb.named('select_inst')[0][1]['asof'] == '2000-01-01T00:00:00'
    assert GalileoSSI._INSTRUMENT_KERNEL is fake_spicedb.kernel_dict
    assert fake_spice.furnished[0] == 'Galileo/SCLK/mk00062a.tsc'
    assert GalileoSSI._initialized


def test_initialize_rejects_an_unexpected_image_shape(fake_spicedb: FakeSpiceDB) -> None:
    fake_spicedb.kernel_dict['INS'][-77036]['MAX_LINE'] = 400
    with pytest.raises(oops.host.HostError,
                       match=r'unsupported SSI image shape \(800, 400\) '
                             r'in the SSI kernel'):
        GalileoSSI.initialize()


def test_initialize_runs_once(fake_spice: FakeSpice) -> None:
    GalileoSSI.initialize()
    GalileoSSI.initialize()
    assert len(fake_spice.solar_system) == 1


@pytest.mark.usefixtures('initialized_ssi')
def test_full_fov() -> None:
    fov = GalileoSSI._FOVS['FULL']
    assert type(fov) is oops.fov.BarrelFOV
    assert tuple(np.asarray(fov.uv_shape.vals)) == (SSI_PIXELS, SSI_PIXELS)
    center = fov.los_from_uv((400., 400.)).unit()
    assert center.vals == pytest.approx((0., 0., 1.))
    scale = SSI_PIXEL_SIZE / SSI_FOCAL_LENGTH
    assert fov.xy_from_uv((401., 400.)).vals == pytest.approx((scale, 0.), rel=1e-6)


@pytest.mark.usefixtures('initialized_ssi')
def test_full_fov_distortion_coefficients() -> None:
    # The kernel's coefficient applies to distances in pixels; the FOV's, in radians
    fov = GalileoSSI._FOVS['FULL']
    scale = SSI_PIXEL_SIZE / SSI_FOCAL_LENGTH
    expected = [1., 0., SSI_DISTORTION / scale**2]
    assert list(fov.coefft_uv_from_xy) == pytest.approx(expected)


@pytest.mark.usefixtures('initialized_ssi')
@pytest.mark.parametrize('mode', ['HIS', 'AI8'])
def test_summed_fovs(mode: str) -> None:
    fov = GalileoSSI._FOVS[mode]
    assert type(fov) is oops.fov.SubsampledFOV
    assert tuple(np.asarray(fov.uv_shape.vals)) == (400, 400)


@pytest.mark.usefixtures('initialized_ssi')
@pytest.mark.parametrize('mode', ['HMA', 'HCA', 'HIM', 'IM8', 'IM4', 'XCM', 'HCJ', 'HCM',
                                  'NONE'])
def test_full_frame_fovs(mode: str) -> None:
    assert GalileoSSI._FOVS[mode] is GalileoSSI._FOVS['FULL']


@pytest.mark.usefixtures('initialized_ssi')
def test_initialize_defines_hd_151288() -> None:
    star = oops.lightsource.star_lookup('HD 151288')
    assert star.magnitude == pytest.approx(8.11)
    assert star.stellar_type == 'K7.5Ve'


def test_initialize_can_repeat_after_reset() -> None:
    GalileoSSI.initialize()
    GalileoSSI._reset()
    GalileoSSI.initialize()
    assert GalileoSSI._initialized


def test_camera_frame_is_defined_once(fake_spice: FakeSpice) -> None:
    GalileoSSI._define_camera_frame()
    GalileoSSI._define_camera_frame()
    assert fake_spice.type1_frames == [('GLL_SCAN_PLATFORM', 40)]


def test_reset_clears_state() -> None:
    GalileoSSI.initialize()
    GalileoSSI._reset()
    assert GalileoSSI._INSTRUMENT_KERNEL is None
    assert GalileoSSI._FOVS == {}
    assert not GalileoSSI._initialized
    assert not _Galileo.initialized

##########################################################################################
# Snapshot construction
##########################################################################################


def test_snapshot_time_is_centered_on_image_time() -> None:
    obs = snapshot()
    tmid = tdb_from_image_time()
    assert obs.time == pytest.approx((tmid - 0.13125, tmid + 0.13125))
    assert obs.texp == pytest.approx(0.2625)
    assert obs.time_from_sclk is False


def test_snapshot_subfields() -> None:
    obs = snapshot()
    assert obs.host == 'Galileo'
    assert obs.instrument == 'SSI'
    assert obs.filter == 'VIOLET'
    assert obs.label_target == 'GANYMEDE'
    assert obs.observation_id == 'G1GSGLOBAL02'
    assert obs.dict['TARGET_NAME'] == 'GANYMEDE'
    assert obs.filepath == FCPath('GO_0017/X.IMG')
    assert obs.basename == 'X.IMG'


def test_snapshot_spice_subfields(fake_spicedb: FakeSpiceDB) -> None:
    obs = snapshot()
    assert obs.spice_kernels == fake_spicedb.used
    assert obs.spice_frame_name == 'GLL_SCAN_PLATFORM'
    assert obs.spice_frame_id == -77001
    assert obs.spice_to_frame == oops.Matrix3.IDENTITY


def test_snapshot_geometry() -> None:
    obs = snapshot()
    assert obs.path.path_id == 'GLL'
    assert obs.frame.frame_id == 'GLL_SCAN_PLATFORM'
    assert obs.target == 'GANYMEDE'
    assert obs.lightsource == 'SUN'


def test_snapshot_uses_cks_for_the_default_frame(fake_spicedb: FakeSpiceDB) -> None:
    snapshot()
    assert fake_spicedb.named('used_basenames')[0][1]['types'] is None


def test_snapshot_custom_frame_reports_no_cks(fake_spice: FakeSpice,
                                              fake_spicedb: FakeSpiceDB) -> None:
    GalileoSSI.initialize()
    # Unlike the fake scan-platform frame, this one is not an identity matrix
    frame = oops.frame.Cmatrix(oops.Matrix3([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),
                               oops.Frame.J2000, frame_id='CUSTOM')
    obs = make_snapshot(label_dict(), frame=frame)
    assert obs.frame is frame
    assert 'CK' not in fake_spicedb.named('used_basenames')[0][1]['types']
    assert fake_spice.type1_frames == []


def test_snapshot_return_all_planets(fake_spicedb: FakeSpiceDB) -> None:
    GalileoSSI.initialize()
    make_snapshot(label_dict(), return_all_planets=True)
    assert 1 in fake_spicedb.named('used_basenames')[0][1]['bodies']


def test_snapshot_cutout_window_trims_the_fov() -> None:
    obs = snapshot(TELEMETRY_FORMAT_ID='IM8')
    assert type(obs.fov) is oops.fov.SliceFOV
    assert tuple(np.asarray(obs.fov.uv_shape.vals)) == (392, 399)
    assert tuple(np.asarray(obs.fov.uv_los.vals)) == pytest.approx((400., 399.))


def test_snapshot_without_cutout_window() -> None:
    obs = snapshot(CUT_OUT_WINDOW=None)
    assert obs.fov is GalileoSSI._FOVS['HIS']


def test_snapshot_supplement_placeholder_window() -> None:
    obs = snapshot(CUT_OUT_WINDOW=np.array([-1, -1, -1, -1]),
                   CUT_OUT_WINDOW_mask=np.array([True, True, True, True]))
    assert obs.fov is GalileoSSI._FOVS['HIS']


def test_snapshot_full_fov_ignores_the_window() -> None:
    GalileoSSI.initialize()
    obs = make_snapshot(label_dict(), full_fov=True)
    assert obs.fov is GalileoSSI._FOVS['HIS']


def test_snapshot_telemetry_format_none() -> None:
    obs = snapshot(TELEMETRY_FORMAT_ID='NONE', CUT_OUT_WINDOW=None)
    assert obs.fov is GalileoSSI._FOVS['NONE']


def test_snapshot_default_telemetry_format() -> None:
    GalileoSSI.initialize()
    obs = make_snapshot(label_dict(TELEMETRY_FORMAT_ID=None, CUT_OUT_WINDOW=None),
                        default_mode='NONE')
    assert obs.fov is GalileoSSI._FOVS['NONE']


@pytest.mark.parametrize('placeholder', ['UNK', 'UNKNOWN', '---', ''])
def test_snapshot_index_null_is_no_value(placeholder: str) -> None:
    obs = snapshot(IMAGE_TIME=placeholder)
    assert obs.time_from_sclk is True


def test_snapshot_trims_the_data() -> None:
    GalileoSSI.initialize()
    obs = make_snapshot(label_dict(), data=image_array())
    assert obs.data.shape == (399, 392)
    assert obs.data[0, 0] == image_array()[1, 0]


@pytest.mark.parametrize(('label_target', 'target'), [
    ('GANYMEDE', 'GANYMEDE'),
    ('BLACK_SKY', 'NONE'),
    ('BLACK SKY', 'NONE'),
    ('NON_SCIENCE', 'NONE'),
    ('NON SCIENCE', 'NONE'),
    ('STAR', 'NONE'),
    ('J RINGS', 'JUPITER_RING_PLANE'),
])
def test_snapshot_target(label_target: str, target: str) -> None:
    obs = snapshot(TARGET_NAME=label_target)
    assert obs.target == target
    assert obs.label_target == label_target
    assert obs.lightsource == 'SUN'


@pytest.mark.parametrize(('label_target', 'star'), [
    ('PLEIADES', 'PLEIADES'),
    ('SIRIUS', 'SIRIUS'),
    ('BET_ARI', 'BET_ARI'),
    ('HD151288', 'HD 151288'),
])
def test_snapshot_star_target(label_target: str, star: str) -> None:
    obs = snapshot(TARGET_NAME=label_target)
    assert obs.target == 'NONE'
    assert obs.lightsource is oops.lightsource.star_lookup(star)


def test_snapshot_overrides(ganymede_path: oops.Path) -> None:
    GalileoSSI.initialize()
    fov = oops.fov.FlatFOV((1.e-5, 1.e-5), (10, 10))
    obs = make_snapshot(label_dict(), target='EUROPA', lightsource='SUN', fov=fov,
                        calibrations=['cal'])
    assert obs.target == 'EUROPA'
    assert obs.fov is fov
    assert obs.calibrations == ['cal']


def test_snapshot_cadence_override() -> None:
    cadence = oops.cadence.SnapCadence(100., 2.)
    obs = make_snapshot(label_dict(), cadence=cadence)
    assert obs.cadence is cadence


def test_snapshot_tracker(ganymede_path: oops.Path) -> None:
    GalileoSSI.initialize()
    obs = make_snapshot(label_dict(), tracker='midtime')
    assert type(obs.frame) is oops.frame.TrackerFrame


def test_snapshot_navigation() -> None:
    GalileoSSI.initialize()
    obs = make_snapshot(label_dict(), navigation=(0., 0.))
    assert type(obs.frame) is oops.frame.Navigation

##########################################################################################
# Times from the spacecraft clock
##########################################################################################


@pytest.mark.parametrize(('frame_duration', 'delay'), [
    (2.333, 0.225), (8.667, 1.149), (15.166, 1.150), (30.333, 2.817), (60.667, 6.483),
], ids=['2.333', '8.667', '15.167', '30.333', '60.667'])
def test_unknown_image_time_comes_from_sclk(frame_duration: float, delay: float) -> None:
    obs = snapshot(IMAGE_TIME='UNK', FRAME_DURATION=frame_duration)
    assert obs.time_from_sclk is True
    assert obs.time[0] == pytest.approx(fake_scs2e(-77, '3496321:0:0:0') + delay,
                                        rel=0, abs=1.e-6)


# The phase-1 flight software ran until RIM 2492218, and the phase-2 software from 3464059
@pytest.mark.parametrize(('sclk', 'shift'), [('02492218.00', 15.167),
                                             ('02492218:00:0:0', 15.167),
                                             ('03464059.00', 0.)],
                         ids=['phase-1', 'phase-1-spice-form', 'phase-2'])
def test_extended_exposure_start(sclk: str, shift: float) -> None:
    normal = snapshot(IMAGE_TIME='UNK', SPACECRAFT_CLOCK_START_COUNT=sclk)
    extended = snapshot(IMAGE_TIME='UNK', SPACECRAFT_CLOCK_START_COUNT=sclk,
                        EXPOSURE_TYPE='EXTENDED')
    assert normal.time[0] - extended.time[0] == pytest.approx(shift)


def test_missing_exposure_type_is_normal() -> None:
    normal = snapshot(IMAGE_TIME='UNK')
    missing = snapshot(IMAGE_TIME='UNK', EXPOSURE_TYPE=None)
    assert missing.time == normal.time


def test_unknown_frame_duration_is_rejected() -> None:
    with pytest.raises(oops.host.HostFormatError,
                       match=r'unrecognized Galileo SSI FRAME_DURATION 20\.0 for '
                             r'GO_0017/X\.IMG'):
        snapshot(IMAGE_TIME='UNK', FRAME_DURATION=20.)

##########################################################################################
# from_file()
##########################################################################################


def test_from_file_pds3(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path))
    assert type(obs) is oops.observation.Snapshot
    assert obs.data.shape == (399, 392)
    assert np.array_equal(obs.data, image_array()[1:400, 0:392])
    assert obs.basename == 'C0349632100R.IMG'
    assert obs.observation_id == 'G1GSGLOBAL02'


def test_from_file_reads_the_window_in_label_order(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, CUT_OUT_WINDOW=[1, 2, 300, 200])
    obs = GalileoSSI.from_file(path)
    assert obs.data.shape == (300, 200)
    assert obs.data[0, 0] == image_array()[0, 1]


def test_from_file_full_fov(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path), full_fov=True)
    assert obs.data.shape == (SSI_PIXELS, SSI_PIXELS)
    assert obs.fov is GalileoSSI._FOVS['HIS']


# A full read retrieves the image file to read its array, so it gives the local path
@pytest.mark.parametrize('given', ['label', 'image'])
def test_from_file_sets_local_paths(tmp_path: pathlib.Path, given: str) -> None:
    label = write_pds3_image(tmp_path)
    image = tmp_path / 'C0349632100R.IMG'
    obs = GalileoSSI.from_file(label if given == 'label' else image)
    assert obs.abspath == image.resolve()
    assert obs.image_url == image.as_posix()


def test_from_file_pds3_filepath_is_the_image(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path), astrometry=True)
    assert obs.filepath.as_posix() == (tmp_path / 'C0349632100R.IMG').as_posix()


_CLOCK: dict[str, float | str | None] = {'RIM': 3496321, 'MOD91': 37, 'MOD10': 4,
                                         'MOD8': 6}


@pytest.mark.parametrize('key', ['MOD10', 'MOD8'])
def test_from_file_pds3_reads_the_vicar_clock(tmp_path: pathlib.Path, key: str) -> None:
    path = write_pds3_image(tmp_path, header=_CLOCK,
                            SPACECRAFT_CLOCK_START_COUNT='03496321.37')
    assert GalileoSSI.from_file(path).dict[key] == _CLOCK[key]


@pytest.mark.parametrize('key', ['MOD10', 'MOD8'])
def test_from_file_vicar_reads_the_clock(tmp_path: pathlib.Path, key: str) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, **_CLOCK), astrometry=True)
    assert obs.dict[key] == _CLOCK[key]


# RIM and MOD91 are wholly inside SPACECRAFT_CLOCK_START_COUNT; MOD10 and MOD8 are not
def test_from_file_vicar_gives_rim_and_mod91_only_in_the_clock_count(
        tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, **_CLOCK), astrometry=True)
    assert {'RIM', 'MOD91'} & obs.dict.keys() == set()
    assert obs.dict['SPACECRAFT_CLOCK_START_COUNT'] == '03496321.37'


def test_from_file_vicar_clock_count_subfield(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, **_CLOCK), astrometry=True)
    assert obs.spacecraft_clock_start_count == '03496321.37'


def test_from_file_pds3_astrometry(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path)
    (tmp_path / 'C0349632100R.IMG').unlink()
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert not hasattr(obs, 'data')
    assert obs.filter == 'VIOLET'
    assert obs.basename == 'C0349632100R.IMG'
    assert 'RIM' not in obs.dict


def test_from_file_vicar(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path))
    assert np.array_equal(obs.data, image_array()[1:400, 0:392])
    assert obs.basename == 'C0349632100R.IMG'
    assert obs.dict['OBSERVATION_ID'] == 'G1GSGLOBAL02'


def test_from_file_vicar_astrometry(tmp_path: pathlib.Path,
                                    monkeypatch: pytest.MonkeyPatch) -> None:
    path = write_vicar_image(tmp_path)

    def fail(*args: Any, **kwargs: Any) -> None:
        raise AssertionError('the VICAR image array was read')

    monkeypatch.setattr(vicar.VicarImage, 'from_file', fail)
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert not hasattr(obs, 'data')
    assert obs.filter == 'VIOLET'


def test_from_file_vicar_image_astrometry(tmp_path: pathlib.Path) -> None:
    path = write_vicar_image(tmp_path)
    vic = cast(vicar.VicarImage, vicar.VicarImage.from_file(path, strict=False))
    vic.filepath = path
    obs = GalileoSSI.from_file(vic, astrometry=True)
    assert not hasattr(obs, 'data')


def test_from_file_vicar_earth_received_time(tmp_path: pathlib.Path) -> None:
    ert = {'ERTYEAR': 1996, 'ERTDAY': 178, 'ERTHOUR': 9, 'ERTMIN': 30, 'ERTSEC': 15,
           'ERTMSEC': 250}
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, **ert), astrometry=True)
    assert obs.dict['EARTH_RECEIVED_TIME'] == '1996-06-26T09:30:15.250'


def test_from_file_vicar_without_earth_received_time(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True)
    assert 'EARTH_RECEIVED_TIME' not in obs.dict


def test_from_file_vicar_label_reads_the_image(tmp_path: pathlib.Path) -> None:
    path = write_vicar_image(tmp_path)
    obs = GalileoSSI.from_file(vicar.VicarLabel(path, strict=False))
    assert np.array_equal(obs.data, image_array()[1:400, 0:392])


def test_from_file_vicar_reads_its_own_header(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True)
    assert obs.filter == 'VIOLET'
    assert obs.observation_id == 'G1GSGLOBAL02'
    assert obs.time[0] == pytest.approx(tdb_from_image_time() - 0.13125, abs=1.e-6)
    assert not obs.time_from_sclk


@pytest.mark.parametrize(('observation_id', 'orbit'), [
    ('G1GSGLOBAL02', 1),
    ('C3CSCRATER01', 3),
    ('E6JSATMOS_01', 6),
    ('20JSSEBTRK12', 20),
])
def test_from_file_vicar_orbit_number(tmp_path: pathlib.Path, observation_id: str,
                                      orbit: int) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, PA=observation_id),
                               astrometry=True)
    assert obs.dict['ORBIT_NUMBER'] == orbit


@pytest.mark.parametrize('observation_id', ['E1LSLUNMAP04____LCLR',
                                            'E2HPTSTARS02*000XCAL', 'OPNAV01A_LCAL',
                                            'GAPSHIPHAS01*016G1IM', '36XSNOTOUR01'])
def test_from_file_vicar_without_orbit_number(tmp_path: pathlib.Path,
                                              observation_id: str) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, PA=observation_id),
                               astrometry=True)
    assert 'ORBIT_NUMBER' not in obs.dict


# Rounding 45.8366 to two places would give 45.84, and truncating it 45.83
@pytest.mark.parametrize('exp', [45.8333, 45.8366])
def test_from_file_vicar_exposure_keeps_its_precision(tmp_path: pathlib.Path,
                                                      exp: float) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, EXP=exp), astrometry=True)
    assert obs.dict['EXPOSURE_DURATION'] == exp


@pytest.mark.parametrize(('rate', 'frame_duration'), [
    (1, 2.333), (2, 8.667), (3, 30.333), (4, 60.667), (5, 15.166),
])
def test_from_file_vicar_frame_duration(tmp_path: pathlib.Path, rate: int,
                                        frame_duration: float) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, RATE=rate), astrometry=True)
    assert obs.dict['FRAME_DURATION'] == frame_duration


@pytest.mark.parametrize(('gain', 'gain_mode'), [
    (1, '400K'), (2, '100K'), (3, '40K'), (4, '10K'),
])
def test_from_file_vicar_gain_mode(tmp_path: pathlib.Path, gain: int,
                                   gain_mode: str) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, GAIN=gain), astrometry=True)
    assert obs.dict['GAIN_MODE_ID'] == gain_mode


def test_from_file_vicar_image_id(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True)
    assert obs.dict['IMAGE_ID'] == 'G1G0002'


def test_from_file_vicar_without_image_id(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, PICNO='?'), astrometry=True)
    assert 'IMAGE_ID' not in obs.dict


@pytest.mark.parametrize(('fibe', 'mofibe', 'exposure_type'), [
    ('1000', None, 'NORMAL'),
    ('1001', None, 'EXTENDED'),
    (None, '001000', 'NORMAL'),
    (None, '001001', 'EXTENDED'),
])
def test_from_file_vicar_exposure_type(tmp_path: pathlib.Path, fibe: str | None,
                                       mofibe: str | None, exposure_type: str) -> None:
    path = write_vicar_image(tmp_path, FIBE=fibe, MOFIBE=mofibe)
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert obs.dict['EXPOSURE_TYPE'] == exposure_type


@pytest.mark.parametrize('missing', [0, -32768])
def test_from_file_vicar_without_image_time(tmp_path: pathlib.Path, missing: int) -> None:
    scet = {f'SCET{field}': missing for field in ('YEAR', 'DAY', 'HOUR', 'MIN', 'SEC',
                                                  'MSEC')}
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path, **scet), astrometry=True)
    assert 'IMAGE_TIME' not in obs.dict
    assert obs.time_from_sclk is True
    assert obs.time[0] == pytest.approx(fake_scs2e(-77, '3496321:0:0:0') + 1.150,
                                        rel=0, abs=1.e-6)


def test_from_file_accepts_a_parsed_label(tmp_path: pathlib.Path) -> None:
    label = Pds3Label(write_pds3_image(tmp_path))
    assert GalileoSSI.from_file(label).data.shape == (399, 392)


def test_from_file_parallel(tmp_path: pathlib.Path) -> None:
    parallel = snapshot()
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True,
                               parallel=parallel)
    assert obs.frame.reference is parallel.frame


def test_from_file_initializes_ssi(tmp_path: pathlib.Path) -> None:
    GalileoSSI.from_file(write_pds3_image(tmp_path), astrometry=True)
    assert GalileoSSI._initialized


def test_host_from_file_dispatches_to_ssi(tmp_path: pathlib.Path) -> None:
    obs = oops.Host.from_file(write_pds3_image(tmp_path), astrometry=True)
    assert obs.instrument == 'SSI'


def test_host_from_file_dispatches_vicar_to_ssi(tmp_path: pathlib.Path) -> None:
    obs = oops.Host.from_file(write_vicar_image(tmp_path), astrometry=True)
    assert obs.instrument == 'SSI'


def test_from_file_rejects_other_hosts(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, SPACECRAFT_NAME='CASSINI ORBITER')
    with pytest.raises(oops.host.WrongHostError, match='not a Galileo SSI PDS3 label'):
        GalileoSSI.from_file(path, astrometry=True)


def test_from_file_rejects_an_unknown_keyword(tmp_path: pathlib.Path) -> None:
    with pytest.raises(TypeError, match="unexpected keyword argument 'astrometric'"):
        GalileoSSI.from_file(write_vicar_image(tmp_path), astrometric=True)


@pytest.mark.parametrize('select', [0, slice(None), 'data'])
def test_from_file_refuses_select(tmp_path: pathlib.Path, select: Any) -> None:
    with pytest.raises(oops.OopsValueError,
                       match='each Galileo SSI image is one observation; select must '
                             'be None'):
        GalileoSSI.from_file(write_vicar_image(tmp_path), select=select)

##########################################################################################
# from_file(): parsed content without a file path
##########################################################################################


def parsed_without_path(tmp_path: pathlib.Path, kind: str) -> Any:
    """The test product parsed from its content, so that it has no file path."""

    label_path = write_pds3_image(tmp_path)
    image_path = tmp_path / 'C0349632100R.IMG'
    if kind == 'label':
        return Pds3Label(label_path.read_text())
    if kind == 'header':
        return vicar.VicarLabel(cast(str, vicar.VicarLabel.read_label(image_path)))
    vic = cast(vicar.VicarImage, vicar.VicarImage.from_file(image_path, strict=False))
    vic.filepath = None
    return vic


@pytest.mark.parametrize(('kind', 'astrometry'), [
    ('label', True), ('header', True), ('image', True), ('image', False),
])
def test_parsed_content_without_a_path_has_no_path_subfields(
        tmp_path: pathlib.Path, kind: str, astrometry: bool) -> None:
    obs = GalileoSSI.from_file(parsed_without_path(tmp_path, kind), astrometry=astrometry)
    assert [name for name in ('filepath', 'basename', 'abspath', 'image_url')
            if hasattr(obs, name)] == []


def test_parsed_image_without_a_path_gives_its_data(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(parsed_without_path(tmp_path, 'image'), full_fov=True)
    assert np.array_equal(obs.data, image_array())


@pytest.mark.parametrize(('kind', 'name'), [('label', 'PDS3 label'),
                                            ('header', 'VICAR header')])
def test_full_read_of_parsed_content_needs_a_path(tmp_path: pathlib.Path, kind: str,
                                                  name: str) -> None:
    with pytest.raises(oops.host.HostFormatError,
                       match=f'no file path is known for the image file of this Galileo '
                             f'SSI {name}'):
        GalileoSSI.from_file(parsed_without_path(tmp_path, kind))

##########################################################################################
# from_file(): what is read
##########################################################################################


def test_image_input_never_reads_the_label(tmp_path: pathlib.Path) -> None:
    write_pds3_image(tmp_path, TARGET_NAME='J RINGS')
    obs = GalileoSSI.from_file(tmp_path / 'C0349632100R.IMG')
    assert obs.dict['TARGET_NAME'] == 'GANYMEDE'


@pytest.mark.parametrize('astrometry', [False, True])
def test_image_input_has_no_label_dict(tmp_path: pathlib.Path, astrometry: bool) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=astrometry)
    assert not hasattr(obs, 'label_dict')


@pytest.mark.parametrize('astrometry', [False, True])
def test_image_input_header_dict(tmp_path: pathlib.Path, astrometry: bool) -> None:
    path = write_vicar_image(tmp_path)
    obs = GalileoSSI.from_file(path, astrometry=astrometry)
    assert obs.header_dict == _ssi_labels.header_json(vicar.VicarLabel(path))


@pytest.mark.parametrize('astrometry', [False, True])
def test_label_input_label_dict(tmp_path: pathlib.Path, astrometry: bool) -> None:
    path = write_pds3_image(tmp_path)
    obs = GalileoSSI.from_file(path, astrometry=astrometry)
    assert obs.label_dict == _ssi_labels.label_json(Pds3Label(path, method='fast')
                                                    .as_dict())


def test_label_input_header_dict(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path))
    image = tmp_path / 'C0349632100R.IMG'
    assert obs.header_dict == _ssi_labels.header_json(vicar.VicarLabel(image))


# rms-vicar's VicarImage replaces the file's HOST with the reading machine's
def test_parsed_image_header_dict_is_the_files(tmp_path: pathlib.Path) -> None:
    path = write_vicar_image(tmp_path)
    vic = cast(vicar.VicarImage, vicar.VicarImage.from_file(path, strict=False))
    vic.filepath = path
    assert GalileoSSI.from_file(vic).header_dict['HOST'] == 'AXP-VMS'


# The processing history repeats TASK, USER, and DAT_TIM for each step
def test_header_dict_keeps_the_whole_processing_history(tmp_path: pathlib.Path) -> None:
    vic = vicar.VicarImage.from_array(image_array())
    vic.label.append([('TASK', 'FIRST'), ('USER', 'U1'), ('TASK', 'SECOND'),
                      ('USER', 'U2')])
    for (key, value) in _VICAR_LABEL.items():
        vic[key] = value
    path = tmp_path / 'C0349632100R.IMG'
    vic.write_file(path)
    header = GalileoSSI.from_file(path, astrometry=True).header_dict
    assert [header.get(key) for key in ('TASK[0]', 'TASK[1]', 'USER[1]')] == [
        'FIRST', 'SECOND', 'U2']


def test_label_is_parsed_without_source_text(tmp_path: pathlib.Path,
                                             monkeypatch: pytest.MonkeyPatch) -> None:
    options: list[dict[str, Any]] = []
    read_fileinfo = GalileoSSI._read_fileinfo

    def spy(fileinfo: Any, **kwargs: Any) -> Any:
        options.append(kwargs['pds3_options'])
        return read_fileinfo(fileinfo, **kwargs)

    monkeypatch.setattr(GalileoSSI, '_read_fileinfo', spy)
    GalileoSSI.from_file(write_pds3_image(tmp_path), astrometry=True,
                         pds3_options={'sources': True})
    assert [option['sources'] for option in options] == [False]


# The parser's own keys, such as "_source", "_fmt", and "objects", are in lower case
def test_label_input_never_gives_the_parsers_keys(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path), astrometry=True,
                               pds3_options={'sources': True})
    keys = {*obs.dict, *obs.label_dict}
    assert [key for key in keys if key != key.upper()] == []


@pytest.mark.parametrize(('input_', 'astrometry'), [('LBL', False), ('LBL', True),
                                                    ('IMG', False), ('IMG', True)],
                         ids=['label', 'label-astrometry', 'image', 'image-astrometry'])
def test_dicts_are_serializable(tmp_path: pathlib.Path, input_: str,
                                astrometry: bool) -> None:
    write_pds3_image(tmp_path, SOURCE_PRODUCT_ID=['A.BSP', 'N/A'])
    obs = GalileoSSI.from_file(tmp_path / f'C0349632100R.{input_}', astrometry=astrometry)
    json.dumps([getattr(obs, name) for name in ('dict', 'label_dict', 'header_dict')
                if hasattr(obs, name)])


@pytest.mark.parametrize('astrometry', [False, True])
def test_label_input_dict_has_no_none(tmp_path: pathlib.Path, astrometry: bool) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path), astrometry=astrometry)
    assert [key for (key, value) in obs.dict.items() if value is None] == []


def test_label_input_dict_has_no_parser_keys(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path), astrometry=True)
    assert [key for key in obs.dict if key != key.upper()] == []


_CONVERTED = ['EXP', 'PA', 'TARGET', 'PICNO', 'FILTER', 'RATE', 'GAIN', 'TLMFMT', 'RIM',
              'MOD91', 'SCETYEAR', 'SCETDAY', 'SCETHOUR', 'SCETMIN', 'SCETSEC',
              'SCETMSEC']


@pytest.mark.parametrize('input_', ['LBL', 'IMG'])
def test_dict_gives_converted_values_only_under_pds3_names(tmp_path: pathlib.Path,
                                                           input_: str) -> None:
    write_pds3_image(tmp_path, header={'TLMFMT': 'HIS'})
    obs = GalileoSSI.from_file(tmp_path / f'C0349632100R.{input_}')
    assert set(_CONVERTED) & obs.dict.keys() == set()


def test_items_have_one_type_from_either_input(tmp_path: pathlib.Path) -> None:
    write_pds3_image(tmp_path)
    from_label = GalileoSSI.from_file(tmp_path / 'C0349632100R.LBL', astrometry=True).dict
    from_image = GalileoSSI.from_file(tmp_path / 'C0349632100R.IMG', astrometry=True).dict
    shared = sorted(from_label.keys() & from_image.keys())
    assert len(shared) >= 10
    assert {key: type(from_label[key]) for key in shared} == {
        key: type(from_image[key]) for key in shared}


def test_label_astrometry_reads_no_header(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path), astrometry=True)
    assert not hasattr(obs, 'header_dict')


def test_label_and_image_give_the_same_data(tmp_path: pathlib.Path) -> None:
    from_label = GalileoSSI.from_file(write_pds3_image(tmp_path), full_fov=True)
    from_image = GalileoSSI.from_file(tmp_path / 'C0349632100R.IMG', full_fov=True)
    assert np.array_equal(from_label.data, from_image.data)


# An archive label can put the image a record early (52 labels of GO_0019) or describe
# the wrong record layout (GO_0018 E4/EUROPA/C0374685452R); the VICAR header is right.
@pytest.mark.parametrize(('record_error', 'samples'), [(-1, SSI_PIXELS), (1, SSI_PIXELS),
                                                       (0, 1000)],
                         ids=['record-early', 'record-late', 'wrong-layout'])
def test_label_input_locates_the_array_by_the_vicar_header(tmp_path: pathlib.Path,
                                                           record_error: int,
                                                           samples: int) -> None:
    path = write_pds3_image(tmp_path, record_error=record_error, samples=samples)
    obs = GalileoSSI.from_file(path, full_fov=True)
    assert np.array_equal(obs.data, image_array())


def test_label_input_finds_the_image_in_another_case(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, image_name='c0349632100r.img')
    obs = GalileoSSI.from_file(path, full_fov=True)
    assert np.array_equal(obs.data, image_array())


# On a case-insensitive file system, the label's own spelling also names the file
def test_label_input_names_the_image_found(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, image_name='c0349632100r.img')
    obs = GalileoSSI.from_file(path)
    assert (tmp_path / obs.basename).samefile(tmp_path / 'c0349632100r.img')


def test_label_input_requires_the_image(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path)
    (tmp_path / 'C0349632100R.IMG').unlink()
    with pytest.raises(FileNotFoundError, match=r'C0349632100R\.IMG'):
        GalileoSSI.from_file(path)


def test_label_input_rejects_an_invalid_image(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path)
    (tmp_path / 'C0349632100R.IMG').write_text('not a VICAR file')
    with pytest.raises(oops.host.VicarHostError, match=r'C0349632100R\.IMG'):
        GalileoSSI.from_file(path)


def test_local_label_astrometry_without_the_image_has_no_abspath(
        tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path)
    (tmp_path / 'C0349632100R.IMG').unlink()
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert (obs.basename, hasattr(obs, 'abspath')) == ('C0349632100R.IMG', False)


_REMOTE = 'https://pds.example.org/volumes/GO_0xxx/GO_0017/G1/GANYMEDE'


def put_in_cache(url: str, content: bytes) -> None:
    """Write `content` where the test's file cache keeps the remote file `url`."""

    local = cast(pathlib.Path, FCPath(url).get_local_path(create_parents=False))
    local.parent.mkdir(parents=True, exist_ok=True)
    local.write_bytes(content)


def cache_remote_product(tmp_path: pathlib.Path, *, image: bool) -> FCPath:
    """Put the test product's label, and its image file if `image` is True, in the local
    cache of the remote directory `_REMOTE`; return the label's URL."""

    write_pds3_image(tmp_path)
    names = ['C0349632100R.LBL'] + (['C0349632100R.IMG'] if image else [])
    for name in names:
        put_in_cache(f'{_REMOTE}/{name}', (tmp_path / name).read_bytes())
    return FCPath(f'{_REMOTE}/C0349632100R.LBL')


# A full read of a remote image file gives the local copy that its retrieval returned
def test_remote_image_full_read_abspath_is_the_retrieved_copy(
        tmp_path: pathlib.Path) -> None:
    cache_remote_product(tmp_path, image=True)
    image = FCPath(f'{_REMOTE}/C0349632100R.IMG')
    obs = GalileoSSI.from_file(image)
    local = cast(pathlib.Path, image.get_local_path(create_parents=False))
    assert obs.abspath == local.resolve()


# An astrometry read does not retrieve the image file, so its Snapshot lacks abspath,
# even when the image file is local or already in the cache
@pytest.mark.parametrize('where', ['local', 'cached'])
def test_label_astrometry_has_no_abspath(tmp_path: pathlib.Path, where: str) -> None:
    label: pathlib.Path | FCPath = write_pds3_image(tmp_path)
    if where == 'cached':
        label = cache_remote_product(tmp_path, image=True)
    obs = GalileoSSI.from_file(label, astrometry=True)
    assert not hasattr(obs, 'abspath')


# The array of a parsed image passed in was read by the caller, not retrieved by this call
def test_parsed_image_has_no_abspath(tmp_path: pathlib.Path) -> None:
    image = vicar.VicarImage(write_vicar_image(tmp_path), strict=False)
    obs = GalileoSSI.from_file(image)
    assert not hasattr(obs, 'abspath')


def test_label_without_an_image_pointer(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path)
    text = path.read_text().replace('^IMAGE ', 'IMAGE_POINTER ')
    path.write_text(text)
    with pytest.raises(oops.host.PdsHostError,
                       match=r'no \^IMAGE pointer in Galileo SSI PDS3 label '
                             r'.*C0349632100R\.LBL'):
        GalileoSSI.from_file(path, astrometry=True)


def test_truncated_image(tmp_path: pathlib.Path) -> None:
    path = write_vicar_image(tmp_path)
    path.write_bytes(path.read_bytes()[:-SSI_PIXELS * 10])
    with pytest.raises(oops.host.VicarHostError, match=r'C0349632100R\.IMG'):
        GalileoSSI.from_file(path)


@pytest.mark.parametrize('parsed', ['label', 'vicar'])
def test_from_file_rejects_another_hosts_parsed_content(tmp_path: pathlib.Path,
                                                        parsed: str) -> None:
    if parsed == 'label':
        content: Any = Pds3Label(write_pds3_image(tmp_path,
                                                  SPACECRAFT_NAME='CASSINI ORBITER'))
    else:
        content = vicar.VicarLabel(write_vicar_image(tmp_path, MISSION='CASSI'))
    with pytest.raises(oops.host.WrongHostError, match='not a Galileo SSI'):
        GalileoSSI.from_file(content, astrometry=True)


def test_from_file_accepts_a_vicar_label(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(vicar.VicarLabel(write_vicar_image(tmp_path)),
                               astrometry=True)
    assert obs.observation_id == 'G1GSGLOBAL02'

def no_pointing(*args: Any, **kwargs: Any) -> None:
    """A frame transform that fails as SPICE does for a time no C kernel covers."""

    raise OSError('SPICE(CKINSUFFDATA) -- sigerr -- insufficient C kernel data')


@pytest.mark.usefixtures('initialized_ssi')
def test_from_file_without_pointing(tmp_path: pathlib.Path,
                                    monkeypatch: pytest.MonkeyPatch) -> None:
    GalileoSSI._define_camera_frame()
    monkeypatch.setattr(oops.Frame.as_frame('GLL_SCAN_PLATFORM'), 'transform_at_time',
                        no_pointing)
    path = write_vicar_image(tmp_path)
    with pytest.raises(oops.OopsSpiceError,
                       match=r'^no Galileo SSI pointing at 1996-06-26T08:46:04\.890 for '
                             r'.*C0349632100R\.IMG: SPICE\(CKINSUFFDATA\)$') as info:
        GalileoSSI.from_file(path, astrometry=True)
    assert info.value.spice_code == 'SPICE(CKINSUFFDATA)'
    # The image file as the filepath subfield names it, with "/" on every platform
    assert info.value.filepath == FCPath(path)
    assert type(info.value.__cause__) is OSError


@pytest.mark.usefixtures('initialized_ssi')
def test_from_file_lets_another_failure_through(tmp_path: pathlib.Path,
                                                monkeypatch: pytest.MonkeyPatch) -> None:
    def fail(*args: Any, **kwargs: Any) -> None:
        raise TypeError('not a SPICE failure')

    GalileoSSI._define_camera_frame()
    monkeypatch.setattr(oops.Frame.as_frame('GLL_SCAN_PLATFORM'), 'transform_at_time',
                        fail)
    with pytest.raises(TypeError, match='not a SPICE failure'):
        GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True)


@pytest.mark.usefixtures('initialized_ssi')
def test_from_file_checks_pointing_in_the_cadence_override(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    GalileoSSI._define_camera_frame()
    times: list[float] = []
    monkeypatch.setattr(oops.Frame.as_frame('GLL_SCAN_PLATFORM'), 'transform_at_time',
                        times.append)
    cadence = oops.cadence.SnapCadence(100., 2.)
    GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True, cadence=cadence)
    assert times == [101.]


# A tracker evaluates the frame once, at the start or end of the exposure that it selects,
# so the pointing is checked there first
@pytest.mark.usefixtures('initialized_ssi', 'ganymede_path')
@pytest.mark.parametrize(('tracker', 'time'), [('start', 100.), ('end', 102.)])
def test_from_file_checks_pointing_at_the_tracker_time(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch, tracker: str,
        time: float) -> None:
    def record(time: Any, quick: Any = None) -> Any:
        times.append(float(oops.Scalar.as_scalar(time).vals))
        return transform(time, quick=quick)

    GalileoSSI._define_camera_frame()
    frame = oops.Frame.as_frame('GLL_SCAN_PLATFORM')
    transform = frame.transform_at_time
    times: list[float] = []
    monkeypatch.setattr(frame, 'transform_at_time', record)
    cadence = oops.cadence.SnapCadence(100., 2.)
    GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True, cadence=cadence,
                         tracker=tracker)
    assert times[0] == time


# A timeshift of the cadence ("T") moves the time at which the frame is evaluated, and so
# does one of the frame ("F"), by the same seconds again if both are given; one of the
# path ("P") leaves it. The letters are read in either case.
@pytest.mark.usefixtures('initialized_ssi')
@pytest.mark.parametrize(('letters', 'time'), [('T', 106.), ('t', 106.), ('F', 106.),
                                               ('f', 106.), ('P', 101.), ('TF', 111.)])
def test_from_file_checks_pointing_at_the_shifted_time(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch, letters: str,
        time: float) -> None:
    # A FrameShift ("F") passes its time as a Scalar, with `quick`
    def record(time: Any, quick: Any = None) -> None:
        times.append(float(oops.Scalar.as_scalar(time).vals))

    GalileoSSI._define_camera_frame()
    times: list[float] = []
    monkeypatch.setattr(oops.Frame.as_frame('GLL_SCAN_PLATFORM'), 'transform_at_time',
                        record)
    cadence = oops.cadence.SnapCadence(100., 2.)
    GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True, cadence=cadence,
                         timeshift=(5., letters))
    assert times == [time]


# A timeshift of the frame alone makes the frame be evaluated later than the midtime, so
# pointing that the C kernels hold at the midtime but not at the shifted time is missing
@pytest.mark.usefixtures('initialized_ssi')
def test_from_file_frame_timeshift_needs_pointing_at_the_shifted_time(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # A FrameShift passes its time as a Scalar, with `quick`
    def pointing_until(time: Any, quick: Any = None) -> None:
        if oops.Scalar.as_scalar(time).vals > 103.:
            no_pointing()

    GalileoSSI._define_camera_frame()
    monkeypatch.setattr(oops.Frame.as_frame('GLL_SCAN_PLATFORM'), 'transform_at_time',
                        pointing_until)
    cadence = oops.cadence.SnapCadence(100., 2.)
    # The message gives the observation's midtime, 101 s, without the frame's seconds
    with pytest.raises(oops.OopsSpiceError,
                       match=r'pointing at 2000-01-01T12:00:36\.816 for .*: '
                             r'SPICE\(CKINSUFFDATA\)$'):
        GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True,
                             cadence=cadence, timeshift=(5., 'F'))


@pytest.mark.usefixtures('initialized_ssi')
def test_from_file_with_a_frame_needs_no_pointing(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    GalileoSSI._define_camera_frame()
    monkeypatch.setattr(oops.Frame.as_frame('GLL_SCAN_PLATFORM'), 'transform_at_time',
                        no_pointing)
    frame = oops.frame.Cmatrix(oops.Matrix3.IDENTITY, oops.Frame.J2000,
                               frame_id='TEST_CMATRIX')
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True, frame=frame)
    assert obs.frame == frame.wayframe


@pytest.mark.usefixtures('initialized_ssi')
def test_from_index_needs_no_pointing(tmp_path: pathlib.Path,
                                      monkeypatch: pytest.MonkeyPatch) -> None:
    GalileoSSI._define_camera_frame()
    monkeypatch.setattr(oops.Frame.as_frame('GLL_SCAN_PLATFORM'), 'transform_at_time',
                        no_pointing)
    assert len(list(GalileoSSI.from_index(write_volume(tmp_path)))) == 3


@pytest.mark.usefixtures('initialized_ssi')
def test_from_index_frame_without_pointing_raises_when_used(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    GalileoSSI._define_camera_frame()
    monkeypatch.setattr(oops.Frame.as_frame('GLL_SCAN_PLATFORM'), 'transform_at_time',
                        no_pointing)
    obs = list(GalileoSSI.from_index(write_volume(tmp_path)))[1]
    with pytest.raises(oops.OopsSpiceError,
                       match=r'^no Galileo SSI pointing at 1996-06-26T08:46:04\.890 for '
                             r'.*C0349632100R\.IMG: SPICE\(CKINSUFFDATA\)$'):
        obs.frame.transform_at_time(obs.midtime)


def test_from_index_frame_with_pointing_is_the_scan_platform(
        tmp_path: pathlib.Path) -> None:
    obs = list(GalileoSSI.from_index(write_volume(tmp_path)))[1]
    assert obs.frame.frame_id == 'GLL_SCAN_PLATFORM'

##########################################################################################
# from_file(): parsing the PDS3 label
##########################################################################################

# The fast parser rejects a quoted value that starts on the line after its keyword; the
# strict parser also rejects an unquoted N/A.
_NEXT_LINE_VALUE = ('NOTE =', '  "FIRST LINE', '  SECOND LINE"')


@pytest.mark.parametrize('extra_lines', [_NEXT_LINE_VALUE,
                                         (*_NEXT_LINE_VALUE, 'SLOPE_FILE_NAME = N/A')],
                         ids=['strict', 'loose'])
def test_from_file_falls_back_to_another_parser(tmp_path: pathlib.Path,
                                                extra_lines: tuple[str, ...]) -> None:
    path = write_pds3_image(tmp_path, extra_lines=extra_lines)
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert 'NOTE' in obs.label_dict


# Only the fast parser accepts a malformed real number
def test_from_file_starts_with_the_requested_parser(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, extra_lines=('NOTE = 1.2.3',))
    with pytest.raises(oops.host.PdsHostError):
        GalileoSSI.from_file(path, astrometry=True, pds3_options={'method': 'strict'})


def test_from_file_passes_another_parser_on(tmp_path: pathlib.Path) -> None:
    with pytest.raises(ValueError, match='bogus'):
        GalileoSSI.from_file(write_pds3_image(tmp_path), astrometry=True,
                             pds3_options={'method': 'bogus'})


def test_from_file_rejects_a_label_no_parser_reads(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, extra_lines=('OBJECT = TABLE',))
    with pytest.raises(oops.host.PdsHostError,
                       match=r'^cannot parse the PDS3 label .*C0349632100R\.LBL: method '
                             r'fast: .*; method strict: .*; method loose: '):
        GalileoSSI.from_file(path, astrometry=True)


# The parsers reject an hour of 25 as a value, not as syntax
def test_from_file_rejects_an_invalid_label_time(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, IMAGE_TIME='1996-06-26T25:46:04.890Z')
    with pytest.raises(oops.host.PdsHostError, match=r'C0349632100R\.LBL'):
        GalileoSSI.from_file(path, astrometry=True)

##########################################################################################
# from_file(): merging the two labels
##########################################################################################


def test_label_only_value_is_merged(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path, PRODUCT_TYPE='REDR'))
    assert obs.dict['PRODUCT_TYPE'] == 'REDR'


def test_header_only_value_is_merged(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path, header={'BOOM': 'N'}))
    assert obs.dict['BOOM'] == 'N'


def test_agreeing_value_is_merged(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path))
    assert obs.dict['OBSERVATION_ID'] == 'G1GSGLOBAL02'


def test_label_wins_a_disagreement(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, header={'TARGET': 'JUPITER'}, TARGET_NAME='J RINGS')
    obs = GalileoSSI.from_file(path)
    assert obs.dict['TARGET_NAME'] == 'J RINGS'


def test_merged_value_sets_the_geometry(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, header={'TARGET': 'JUPITER'}, TARGET_NAME='J RINGS')
    obs = GalileoSSI.from_file(path)
    assert obs.target == 'JUPITER_RING_PLANE'


def test_header_exposure_wins(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, header={'EXP': 45.8333}, EXPOSURE_DURATION=45.83)
    obs = GalileoSSI.from_file(path)
    assert obs.texp == 0.0458333


# The rule names a source, not the value with more digits
def test_header_distance_wins_with_fewer_digits(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, header={'TARGET_CENTER_DISTANCE': 665796.},
                            TARGET_CENTER_DISTANCE=665795.8)
    obs = GalileoSSI.from_file(path)
    assert obs.dict['TARGET_CENTER_DISTANCE'] == 665796.


def test_label_exposure_fills_a_header_without_one(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, header={'EXP': None}, EXPOSURE_DURATION=45.83)
    assert GalileoSSI.from_file(path).texp == pytest.approx(0.04583, rel=1.e-12)


def test_label_number_wins_a_disagreement(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, header={'PLRANGE': 1689750.},
                            CENTRAL_BODY_DISTANCE=1689682.)
    obs = GalileoSSI.from_file(path)
    assert obs.dict['CENTRAL_BODY_DISTANCE'] == 1689682.


def test_label_orbit_number_outside_the_tour_is_dropped(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, header={'PA': 'E1WNMESON-01WDRK'},
                            OBSERVATION_ID='E1WNMESON-01WDRK', ORBIT_NUMBER=12613)
    obs = GalileoSSI.from_file(path)
    assert 'ORBIT_NUMBER' not in obs.dict


def test_label_dict_keeps_a_dropped_value(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, header={'PA': 'E1WNMESON-01WDRK'},
                            OBSERVATION_ID='E1WNMESON-01WDRK', ORBIT_NUMBER=12613)
    obs = GalileoSSI.from_file(path)
    assert obs.label_dict['ORBIT_NUMBER'] == 12613

##########################################################################################
# from_file(): missing values
##########################################################################################


def test_label_null_is_not_merged(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path, SMEAR_AZIMUTH='UNK'))
    assert 'SMEAR_AZIMUTH' not in obs.dict


def test_label_null_stays_in_the_label_dict(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path, SMEAR_AZIMUTH='UNK'))
    assert obs.label_dict['SMEAR_AZIMUTH'] == 'UNK'


@pytest.mark.parametrize('key', ['SMRAZ', 'SMEAR_AZIMUTH'])
def test_header_fill_is_not_merged(tmp_path: pathlib.Path, key: str) -> None:
    path = write_vicar_image(tmp_path, SMRAZ=-999.)
    assert key not in GalileoSSI.from_file(path, astrometry=True).dict


def test_header_fill_stays_in_the_header_dict(tmp_path: pathlib.Path) -> None:
    path = write_vicar_image(tmp_path, SMRAZ=-999.)
    assert GalileoSSI.from_file(path, astrometry=True).header_dict['SMRAZ'] == -999.


def test_header_input_gives_no_placeholders(tmp_path: pathlib.Path) -> None:
    scet = {f'SCET{field}': -32768 for field in ('YEAR', 'DAY', 'HOUR', 'MIN', 'SEC',
                                                 'MSEC')}
    path = write_vicar_image(tmp_path, PICNO='?', PA='OPNAV01A_LCAL', SMRAZ=-999., **scet)
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert [key for (key, value) in obs.dict.items() if value in (None, 'UNK')] == []


@pytest.mark.parametrize(('changes', 'message'), [
    ({'EXPOSURE_DURATION': None}, 'exposure unknown: no EXPOSURE_DURATION'),
    ({'IMAGE_TIME': 'UNK', 'SPACECRAFT_CLOCK_START_COUNT': None},
     'time unknown: no SPACECRAFT_CLOCK_START_COUNT'),
    ({'IMAGE_TIME': 'UNK', 'FRAME_DURATION': None}, 'time unknown: no FRAME_DURATION'),
    ({'TELEMETRY_FORMAT_ID': None}, 'telemetry format unknown: no TELEMETRY_FORMAT_ID'),
], ids=['exposure', 'clock-count', 'frame-duration', 'telemetry-format'])
def test_from_file_requires_a_value(tmp_path: pathlib.Path, changes: dict[str, Any],
                                    message: str) -> None:
    path = write_pds3_image(tmp_path, **changes)
    with pytest.raises(oops.host.HostFormatError,
                       match=f'Galileo SSI {message} for .*C0349632100R\\.IMG'):
        GalileoSSI.from_file(path, astrometry=True)


def test_from_file_rejects_an_unknown_telemetry_format(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path, TELEMETRY_FORMAT_ID='XYZ')
    with pytest.raises(oops.host.HostFormatError,
                       match="unrecognized Galileo SSI TELEMETRY_FORMAT_ID 'XYZ'"):
        GalileoSSI.from_file(path, astrometry=True)


@pytest.mark.parametrize('name', ['filter', 'label_target', 'observation_id',
                                  'spacecraft_clock_start_count'])
def test_missing_value_has_no_subfield(tmp_path: pathlib.Path, name: str) -> None:
    path = write_pds3_image(tmp_path, FILTER_NAME=None, TARGET_NAME=None,
                            OBSERVATION_ID=None, SPACECRAFT_CLOCK_START_COUNT='UNK')
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert not hasattr(obs, name)


def test_label_clock_count_subfield(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path), astrometry=True)
    assert obs.spacecraft_clock_start_count == SCLK_COUNT

##########################################################################################
# from_index()
##########################################################################################

_INDEX_COLUMNS = [
    ('VOLUME_ID', 'CHARACTER', 7),
    ('FILE_SPECIFICATION_NAME', 'CHARACTER', 30),
    ('SPACECRAFT_CLOCK_START_COUNT', 'CHARACTER', 11),
    ('MISSION_NAME', 'CHARACTER', 7),
    ('INSTRUMENT_ID', 'CHARACTER', 3),
    ('OBSERVATION_ID', 'CHARACTER', 16),
    ('TARGET_NAME', 'CHARACTER', 10),
    ('IMAGE_TIME', 'CHARACTER', 24),
    ('FILTER_NAME', 'CHARACTER', 7),
    ('EXPOSURE_DURATION', 'ASCII_REAL', 9),
    ('FRAME_DURATION', 'ASCII_REAL', 9),
]

_INDEX_ROWS = [
    ('GO_0017', 'J0/OPNAV/C0346405900R.LBL', '03464059.00', 'GALILEO', 'SSI',
     'J0GSOPNAV115', 'GANYMEDE', '1996-06-03T17:05:38.015Z', 'CLEAR', 12.5, 8.667),
    ('GO_0017', 'G1/GANYMEDE/C0349632100R.LBL', SCLK_COUNT, 'GALILEO', 'SSI',
     'G1GSGLOBAL02', 'GANYMEDE', IMAGE_TIME, 'VIOLET', EXPOSURE_MS, 15.166),
    ('GO_0017', 'G1/GANYMEDE/C0349632200R.LBL', '03496322.00', 'GALILEO', 'SSI',
     'G1GSGLOBAL02', 'GANYMEDE', '1996-06-26T08:47:05.557Z', 'CLEAR', 0., 15.166),
]

# The supplement declares a CUT_OUT_WINDOW of -1, -1, -1, -1 as not applicable
_SUPPLEMENT_COLUMNS = [
    ('VOLUME_ID', 'CHARACTER', 7),
    ('FILE_SPECIFICATION_NAME', 'CHARACTER', 30),
    ('TELEMETRY_FORMAT_ID', 'CHARACTER', 3),
    ('CUT_OUT_WINDOW', 'ASCII_INTEGER', 4, -1),
    ('EXPOSURE_TYPE', 'CHARACTER', 8),
    ('SPACECRAFT_CLOCK_START_COUNT', 'CHARACTER', 15),
]

_SUPPLEMENT_ROWS = [
    ('GO_0017', 'J0/OPNAV/C0346405900R.LBL', 'IM8', [-1, -1, -1, -1], 'NORMAL',
     '03464059:00:0:0'),
    ('GO_0017', 'G1/GANYMEDE/C0349632100R.LBL', 'HIS', [2, 1, 399, 392], 'NORMAL',
     '03496321:00:0:0'),
    ('GO_0017', 'G1/GANYMEDE/C0349632200R.LBL', 'HIS', [2, 1, 399, 392], 'NORMAL',
     '03496322:00:0:0'),
]


def write_table(directory: pathlib.Path, name: str, object_name: str,
                columns: list[tuple[Any, ...]], rows: list[tuple[Any, ...]],
                **keywords: Any) -> pathlib.Path:
    """Write an ASCII PDS3 table and its label; return the label's path.

    Each column is (name, data type, width), or (name, data type, width, constant) for a
    column that declares a NOT_APPLICABLE_CONSTANT. An ASCII_INTEGER column holds a list
    of four values, each in a field of `width` bytes; an ASCII_REAL column holds one
    value.
    """

    def fields(value: Any, kind: str, width: int) -> list[str]:
        if kind == 'CHARACTER':
            return ['"' + str(value).ljust(width) + '"']
        if kind == 'ASCII_INTEGER':
            return [f'{item:{width}d}' for item in value]
        return [f'{value:{width}.3f}']

    lines = [','.join(item for (value, (_, kind, width, *_)) in zip(row, columns,
                                                                    strict=True)
                      for item in fields(value, kind, width)) + '\r\n' for row in rows]
    (directory / f'{name}.tab').write_text(''.join(lines), newline='')
    row_bytes = len(lines[0])

    text = ['PDS_VERSION_ID = PDS3', 'RECORD_TYPE = FIXED_LENGTH',
            f'RECORD_BYTES = {row_bytes}', f'FILE_RECORDS = {len(rows)}',
            f'^{object_name} = "{name}.tab"']
    text += pds3_lines(keywords)
    text += [f'OBJECT = {object_name}', '  INTERCHANGE_FORMAT = ASCII',
             f'  ROWS = {len(rows)}', f'  COLUMNS = {len(columns)}',
             f'  ROW_BYTES = {row_bytes}']
    start = 1
    for (column, kind, width, *constant) in columns:
        text += ['  OBJECT = COLUMN', f'    NAME = {column}', f'    DATA_TYPE = {kind}']
        text += [f'    NOT_APPLICABLE_CONSTANT = {value}' for value in constant]
        if kind == 'CHARACTER':
            text += [f'    START_BYTE = {start + 1}', f'    BYTES = {width}']
            start += width + 3
        elif kind == 'ASCII_INTEGER':
            text += [f'    START_BYTE = {start}', f'    BYTES = {4 * width + 3}',
                     '    ITEMS = 4', f'    ITEM_BYTES = {width}',
                     f'    ITEM_OFFSET = {width + 1}']
            start += 4 * (width + 1)
        else:
            text += [f'    START_BYTE = {start}', f'    BYTES = {width}']
            start += width + 1
        text += ['  END_OBJECT = COLUMN']
    text += [f'END_OBJECT = {object_name}']
    return write_text(directory / f'{name}.lbl', text)


_INDEX_KEYWORDS = {'MISSION_NAME': 'GALILEO', 'INSTRUMENT_HOST_NAME': 'GALILEO',
                   'INSTRUMENT_ID': 'SSI'}


def write_index(directory: pathlib.Path,
                rows: list[tuple[Any, ...]] = _INDEX_ROWS) -> pathlib.Path:
    """Write an SSI index, by default of three rows, and its label; return the label's
    path."""

    return write_table(directory, 'GO_0017_index', 'IMAGE_INDEX_TABLE', _INDEX_COLUMNS,
                       rows, **_INDEX_KEYWORDS)


def write_supplement(directory: pathlib.Path,
                     rows: list[tuple[Any, ...]] = _SUPPLEMENT_ROWS) -> pathlib.Path:
    """Write the supplemental index for :func:`write_index`; return its label's path."""

    return write_table(directory, 'GO_0017_supplemental_index',
                       'SUPPLEMENTAL_INDEX_TABLE', _SUPPLEMENT_COLUMNS, rows,
                       **_INDEX_KEYWORDS)


def write_volume(directory: pathlib.Path) -> pathlib.Path:
    """Write the SSI index and its supplemental index; return the index label's path."""

    write_supplement(directory)
    return write_index(directory)


def index_rows(index: pathlib.Path) -> list[dict[str, Any]]:
    """The rows of an SSI index alone, as dictionaries, for `row_dicts`."""

    return list(oops.Host._index_rows(index, supplement=None))


def test_index_and_supplement_are_paired_by_position() -> None:
    assert GalileoSSI.MERGE_SUPPLEMENT_BY_FILE is False


# The pairing compares the files of the two rows, so a supplement without them is refused
# at once rather than paired unchecked
def test_from_index_needs_the_file_specification_in_the_supplement(
        tmp_path: pathlib.Path) -> None:
    columns = [_SUPPLEMENT_COLUMNS[0], *_SUPPLEMENT_COLUMNS[2:]]
    rows = [(row[0], *row[2:]) for row in _SUPPLEMENT_ROWS]
    write_table(tmp_path, 'GO_0017_supplemental_index', 'SUPPLEMENTAL_INDEX_TABLE',
                columns, rows, **_INDEX_KEYWORDS)
    with pytest.raises(oops.host.HostIndexError, match='FILE_SPECIFICATION_NAME'):
        GalileoSSI.from_index(write_index(tmp_path))


def long_volume(directory: pathlib.Path, count: int) -> pathlib.Path:
    """Write an index and supplement of `count` rows, each a copy of the second row of
    :data:`_INDEX_ROWS` for another image; return the index label's path."""

    def renamed(row: tuple[Any, ...], number: int) -> tuple[Any, ...]:
        return (row[0], f'G1/GANYMEDE/C03496{number:04d}R.LBL', *row[2:])

    write_supplement(directory, [renamed(_SUPPLEMENT_ROWS[1], k) for k in range(count)])
    return write_index(directory, [renamed(_INDEX_ROWS[1], k) for k in range(count)])


def test_from_index_reads_no_table_at_call_time(tmp_path: pathlib.Path) -> None:
    """The call reads only the labels: with both table files gone it still returns, and
    the first step of the iteration finds them missing."""

    index = write_volume(tmp_path)
    for table in tmp_path.glob('*.tab'):
        table.unlink()
    snapshots = GalileoSSI.from_index(index)
    with pytest.raises(FileNotFoundError):
        next(snapshots)


def test_from_index_yields_rows_before_a_later_missing_row(
        tmp_path: pathlib.Path) -> None:
    """The rows are read as the iteration proceeds: a row that the label promises but
    the table file lacks, beyond the first chunks of rows read, is found only after the
    rows before it are yielded."""

    count = 2 * CHUNK_ROWS + 1
    index = long_volume(tmp_path, count)
    table = tmp_path / 'GO_0017_index.tab'
    table.write_bytes(table.read_bytes().rpartition(b'"GO_0017"')[0])
    snapshots = GalileoSSI.from_index(index)
    assert next(snapshots).basename == 'C034960000R.IMG'
    assert len(list(itertools.islice(snapshots, count - 2))) == count - 2
    with pytest.raises(oops.host.HostIndexError,
                       match=r'invalid index table .*GO_0017_index\.lbl'):
        next(snapshots)


def test_from_index_builds_each_row_when_reached(tmp_path: pathlib.Path,
                                                 monkeypatch: pytest.MonkeyPatch) -> None:
    """The result is an iterator that builds each row's Snapshot when it reaches the row
    and yields the rows in index order."""

    built: list[Any] = []
    make_snapshot = GalileoSSI._make_snapshot

    def spy(dict_: dict[str, Any], **kwargs: Any) -> Any:
        built.append(kwargs['filepath'])
        return make_snapshot(dict_, **kwargs)

    monkeypatch.setattr(GalileoSSI, '_make_snapshot', spy)
    snapshots = GalileoSSI.from_index(write_volume(tmp_path))
    assert isinstance(snapshots, Iterator)
    assert built == []
    first = next(snapshots)
    assert len(built) == 1
    assert [first.basename] + [obs.basename for obs in snapshots] == [
        'C0346405900R.IMG', 'C0349632100R.IMG', 'C0349632200R.IMG']


def test_from_index_raises_a_row_error_when_reached(tmp_path: pathlib.Path) -> None:
    """A row that cannot give a Snapshot raises when the iteration reaches it, after the
    rows before it have been yielded."""

    index = write_index(tmp_path)
    rows = index_rows(index)
    del rows[1]['EXPOSURE_DURATION']
    snapshots = GalileoSSI.from_index(index, row_dicts=rows)
    assert next(snapshots).basename == 'C0346405900R.IMG'
    with pytest.raises(oops.host.HostFormatError,
                       match=r'exposure unknown: no EXPOSURE_DURATION for '
                             r'.*C0349632100R\.IMG'):
        next(snapshots)


def test_from_index_clock_error_names_the_image(tmp_path: pathlib.Path,
                                                monkeypatch: pytest.MonkeyPatch) -> None:
    """A row timed from a clock count that the clock kernel cannot convert raises an
    error naming that row's image file."""

    def scs2e(spacecraft: int, count: str) -> float:
        raise RuntimeError('SPICE(VALUEOUTOFRANGE) -- the count is outside the clock')

    monkeypatch.setattr(cspyce, 'scs2e', scs2e)
    index = write_index(tmp_path)
    rows = index_rows(index)[1:2]
    rows[0]['IMAGE_TIME'] = 'UNK'
    snapshots = GalileoSSI.from_index(index, row_dicts=rows)
    with pytest.raises(oops.OopsSpiceError,
                       match=r"^no time for Galileo clock count '03496321\.00' for "
                             r".*C0349632100R\.IMG: SPICE\(VALUEOUTOFRANGE\)$"):
        next(snapshots)


def test_from_index_keeps_zero_exposures(tmp_path: pathlib.Path) -> None:
    snapshots = list(GalileoSSI.from_index(write_volume(tmp_path)))
    assert [obs.basename for obs in snapshots] == ['C0346405900R.IMG', 'C0349632100R.IMG',
                                                   'C0349632200R.IMG']
    assert snapshots[2].texp == 0.


def test_from_index_rows(tmp_path: pathlib.Path) -> None:
    snapshots = list(GalileoSSI.from_index(write_volume(tmp_path)))
    assert [str(obs.filepath) for obs in snapshots] == [
        (tmp_path / 'GO_0017' / spec).as_posix()
        for spec in ('J0/OPNAV/C0346405900R.IMG', 'G1/GANYMEDE/C0349632100R.IMG',
                     'G1/GANYMEDE/C0349632200R.IMG')]
    assert [obs.filter for obs in snapshots] == ['CLEAR', 'VIOLET', 'CLEAR']
    assert [obs.observation_id for obs in snapshots] == ['J0GSOPNAV115', 'G1GSGLOBAL02',
                                                         'G1GSGLOBAL02']
    assert snapshots[1].time[0] == pytest.approx(tdb_from_image_time() - 0.13125)


def test_from_index_names_the_image_in_full(tmp_path: pathlib.Path) -> None:
    image = tmp_path / 'GO_0017/G1/GANYMEDE/C0349632100R.IMG'
    obs = list(GalileoSSI.from_index(write_volume(tmp_path)))[1]
    assert (obs.basename, obs.image_url) == ('C0349632100R.IMG', image.as_posix())


# An index read never retrieves an image file, so no row has an abspath, even when its
# image file is local or already in the cache
@pytest.mark.parametrize('where', ['local', 'cached'])
def test_from_index_has_no_abspath(tmp_path: pathlib.Path, where: str) -> None:
    local = write_volume(tmp_path)
    index: pathlib.Path | str = local
    if where == 'local':
        image = tmp_path / 'GO_0017/G1/GANYMEDE/C0349632100R.IMG'
        image.parent.mkdir(parents=True)
        image.write_bytes(b'')
    else:
        put_in_cache(f'{_REMOTE}/C0349632100R.IMG', b'')
        index = 'https://pds.example.org/metadata/GO_0xxx/GO_0017/GO_0017_index.lbl'
    snapshots = GalileoSSI.from_index(index, row_dicts=index_rows(local),
                                      label=Pds3Label(local))
    assert [hasattr(obs, 'abspath') for obs in snapshots] == [False, False, False]


# An index whose volumes are absent leaves the file system as it was
@pytest.mark.parametrize('remote', [False, True], ids=['local', 'remote'])
def test_from_index_creates_no_directories(tmp_path: pathlib.Path, remote: bool) -> None:
    (tmp_path / 'metadata').mkdir()
    local = write_volume(tmp_path / 'metadata')
    index: pathlib.Path | FCPath = local
    rows = index_rows(local)
    if remote:
        index = FCPath('https://pds.example.org/metadata/GO_0017_index.lbl')
    before = sorted(tmp_path.rglob('*'))
    snapshots = list(GalileoSSI.from_index(index, row_dicts=rows,
                                           label=Pds3Label(local)))
    assert [hasattr(obs, 'abspath') for obs in snapshots] == [False, False, False]
    assert sorted(tmp_path.rglob('*')) == before


# The image paths of index rows come from the location of the index
def test_from_index_finds_the_volumes_beside_the_metadata(tmp_path: pathlib.Path) -> None:
    metadata = tmp_path / 'metadata' / 'GO_0xxx' / 'GO_0017'
    metadata.mkdir(parents=True)
    obs = next(GalileoSSI.from_index(write_volume(metadata)))
    assert obs.filepath == FCPath(tmp_path / 'volumes/GO_0xxx/GO_0017/J0/OPNAV/'
                                  'C0346405900R.IMG')


def test_from_index_parses_the_index_label_once(tmp_path: pathlib.Path,
                                                monkeypatch: pytest.MonkeyPatch) -> None:
    index = write_volume(tmp_path)
    parsed: list[FCPath] = []
    parse = oops.Host._parse_index_label

    def spy(filepath: Any, *args: Any, **kwargs: Any) -> Any:
        parsed.append(FCPath(filepath))
        return parse(filepath, *args, **kwargs)

    monkeypatch.setattr(oops.Host, '_parse_index_label', spy)
    GalileoSSI.from_index(index)
    assert parsed.count(FCPath(index)) == 1


def test_from_index_requires_the_supplement_beside_it(tmp_path: pathlib.Path) -> None:
    with pytest.raises(FileNotFoundError, match=r'GO_0017_supplemental_index\.lbl'):
        GalileoSSI.from_index(write_index(tmp_path))


def test_from_index_with_supplement(tmp_path: pathlib.Path) -> None:
    snapshots = list(GalileoSSI.from_index(write_index(tmp_path),
                                           supplement=write_supplement(tmp_path)))
    assert snapshots[0].fov is GalileoSSI._FOVS['IM8']
    assert tuple(snapshots[1].fov.uv_shape.vals) == (392, 399)


def test_from_index_keeps_the_index_values(tmp_path: pathlib.Path) -> None:
    snapshots = list(GalileoSSI.from_index(write_index(tmp_path),
                                           supplement=write_supplement(tmp_path)))
    assert snapshots[1].dict['SPACECRAFT_CLOCK_START_COUNT'] == SCLK_COUNT


def test_from_index_finds_the_supplement(tmp_path: pathlib.Path) -> None:
    write_supplement(tmp_path)
    snapshots = list(GalileoSSI.from_index(write_index(tmp_path)))
    assert snapshots[0].fov is GalileoSSI._FOVS['IM8']


def test_from_index_requires_a_named_supplement(tmp_path: pathlib.Path) -> None:
    with pytest.raises(FileNotFoundError, match=r'missing_supplemental_index\.lbl'):
        GalileoSSI.from_index(write_index(tmp_path),
                              supplement=tmp_path / 'missing_supplemental_index.lbl')


def test_from_index_rejects_another_host(tmp_path: pathlib.Path) -> None:
    path = write_table(tmp_path, 'GO_0017_index', 'IMAGE_INDEX_TABLE', _INDEX_COLUMNS,
                       _INDEX_ROWS, MISSION_NAME='CASSINI',
                       INSTRUMENT_HOST_NAME='CASSINI', INSTRUMENT_ID='ISS')
    with pytest.raises(oops.host.WrongHostError, match='not a Galileo SSI index'):
        GalileoSSI.from_index(path)


# The rows are paired by position, so a supplement in another order raises at the first
# row whose files differ, after the rows before it are yielded
def test_from_index_rejects_a_supplement_in_another_order(tmp_path: pathlib.Path) -> None:
    rows = [_SUPPLEMENT_ROWS[0], _SUPPLEMENT_ROWS[2], _SUPPLEMENT_ROWS[1]]
    snapshots = GalileoSSI.from_index(write_index(tmp_path),
                                      supplement=write_supplement(tmp_path, rows))
    assert next(snapshots).basename == 'C0346405900R.IMG'
    with pytest.raises(oops.host.HostIndexError,
                       match=r"index row 1 gives '.*C0349632100R\.LBL' but the "
                             r"supplement row in the same position gives "
                             r"'.*C0349632200R\.LBL'"):
        next(snapshots)


def test_from_index_rejects_a_short_supplement(tmp_path: pathlib.Path) -> None:
    supplement = write_supplement(tmp_path, _SUPPLEMENT_ROWS[:2])
    snapshots = GalileoSSI.from_index(write_index(tmp_path), supplement=supplement)
    assert len(list(itertools.islice(snapshots, 2))) == 2
    with pytest.raises(oops.host.HostIndexError,
                       match=r'the supplement ends before index row 2'):
        next(snapshots)


def test_from_index_rejects_a_long_supplement(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(
        write_index(tmp_path, _INDEX_ROWS[:2]), supplement=write_supplement(tmp_path))
    assert len(list(itertools.islice(snapshots, 2))) == 2
    with pytest.raises(oops.host.HostIndexError,
                       match='the supplement has more rows than the index'):
        next(snapshots)


def test_from_index_filter(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_volume(tmp_path),
                                      filter={'FILTER_NAME': 'VIOLET'})
    assert [obs.filter for obs in snapshots] == ['VIOLET']


def test_from_index_filter_uses_the_supplement(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_volume(tmp_path),
                                      filter={'TELEMETRY_FORMAT_ID': 'IM8'})
    assert [obs.basename for obs in snapshots] == ['C0346405900R.IMG']


def test_from_index_refuses_an_unknown_filter_column(tmp_path: pathlib.Path) -> None:
    """The call itself raises, before any row is read."""

    with pytest.raises(oops.OopsValueError, match='NOT_A_COLUMN'):
        GalileoSSI.from_index(write_volume(tmp_path), filter={'NOT_A_COLUMN': 1})


def test_from_index_uses_given_rows(tmp_path: pathlib.Path) -> None:
    index = write_index(tmp_path)
    snapshots = GalileoSSI.from_index(index, row_dicts=index_rows(index)[:1])
    assert [obs.basename for obs in snapshots] == ['C0346405900R.IMG']


def test_from_index_iterates_given_rows_lazily(tmp_path: pathlib.Path) -> None:
    """Rows given as an iterator are taken one at a time, as the iteration needs them."""

    index = write_index(tmp_path)
    taken: list[str] = []

    def rows() -> Iterator[dict[str, Any]]:
        for row in index_rows(index):
            taken.append(row['FILE_SPECIFICATION_NAME'])
            yield row

    snapshots = GalileoSSI.from_index(index, row_dicts=rows())
    next(snapshots)
    assert taken == ['J0/OPNAV/C0346405900R.LBL']


def test_from_index_names_the_image_in_the_label_case(tmp_path: pathlib.Path) -> None:
    index = write_index(tmp_path)
    rows = index_rows(index)
    rows[0]['FILE_SPECIFICATION_NAME'] = 'j0/opnav/c0346405900r.lbl'
    [obs] = GalileoSSI.from_index(index, row_dicts=rows[:1])
    assert obs.basename == 'c0346405900r.img'


_RATIOS = ('ENCODING_COMPRESSION_RATIO', 'ENCODING_MIN_COMPRESSION_RATIO',
           'ENCODING_MAX_COMPRESSION_RATIO')
_ROTATED = ('ENCODING_MIN_COMPRESSION_RATIO', 'ENCODING_MAX_COMPRESSION_RATIO',
            'ENCODING_COMPRESSION_RATIO')


# The current indexes, whose labels name the ratio columns in the rotated order, give each
# value under the name of the ratio it holds, also when the label is made again with a
# later creation time, as every earlier regeneration was; an index corrected by renaming
# the columns keeps the names it gives
@pytest.mark.parametrize(('order', 'created', 'ratios'), [
    (_ROTATED, '2025-10-30T13:48:44', (1.5, 2.5, 3.5)),
    (_RATIOS, '2025-10-30T13:48:44', (1.5, 2.5, 3.5)),
    (_ROTATED, '2026-12-01T00:00:00', (1.5, 2.5, 3.5)),
], ids=['rotated', 'renamed', 'made-again'])
def test_from_index_ratio_columns(tmp_path: pathlib.Path, order: tuple[str, ...],
                                  created: str, ratios: tuple[float, ...]) -> None:
    columns = _INDEX_COLUMNS + [(name, 'ASCII_REAL', 8) for name in order]
    rows = [(*row, 1.5, 2.5, 3.5) for row in _INDEX_ROWS]
    write_supplement(tmp_path)
    index = write_table(tmp_path, 'GO_0017_index', 'IMAGE_INDEX_TABLE', columns, rows,
                        PRODUCT_CREATION_TIME=created, **_INDEX_KEYWORDS)
    obs = next(GalileoSSI.from_index(index))
    assert tuple(obs.dict[name] for name in _RATIOS) == ratios


def test_from_index_clock_count_subfield(tmp_path: pathlib.Path) -> None:
    snapshots = list(GalileoSSI.from_index(write_volume(tmp_path)))
    assert snapshots[1].spacecraft_clock_start_count == SCLK_COUNT


def test_from_index_fov_override(tmp_path: pathlib.Path) -> None:
    fov = oops.fov.FlatFOV((1.e-5, 1.e-5), (10, 10))
    snapshots = GalileoSSI.from_index(write_volume(tmp_path), fov=fov)
    assert all(obs.fov is fov for obs in snapshots)


def test_from_index_parallel(tmp_path: pathlib.Path) -> None:
    parallel = snapshot()
    snapshots = list(GalileoSSI.from_index(write_volume(tmp_path), parallel=parallel))
    assert all(type(obs.frame) is oops.frame.Cmatrix for obs in snapshots)
    assert all(obs.frame.reference is parallel.frame for obs in snapshots)


def test_host_from_index_dispatches_to_ssi(tmp_path: pathlib.Path) -> None:
    snapshots = oops.Host.from_index(write_volume(tmp_path))
    assert [obs.instrument for obs in snapshots] == ['SSI', 'SSI', 'SSI']


# The supplement is required, so None is refused rather than taken as "search for it"
def test_from_index_refuses_no_supplement(tmp_path: pathlib.Path) -> None:
    index = write_volume(tmp_path)
    with pytest.raises(oops.host.HostIndexError,
                       match=r'^Galileo SSI index needs its supplemental index: '
                             + re.escape(str(FCPath(index))) + '$'):
        GalileoSSI.from_index(index, supplement=None)


def test_from_index_rows_need_no_supplement(tmp_path: pathlib.Path) -> None:
    index = write_index(tmp_path)
    assert len(list(GalileoSSI.from_index(index, row_dicts=index_rows(index),
                                          supplement=None))) == 3


def test_from_index_given_rows_need_the_label(tmp_path: pathlib.Path) -> None:
    """The label decides whether the ratio columns are rotated, so it is read."""

    rows = index_rows(write_index(tmp_path))
    with pytest.raises(FileNotFoundError, match=r'missing\.lbl'):
        GalileoSSI.from_index(tmp_path / 'missing.lbl', row_dicts=rows)


def test_from_index_uses_a_given_label(tmp_path: pathlib.Path,
                                       monkeypatch: pytest.MonkeyPatch) -> None:
    """A given label is used as it is and the index label is not parsed again."""

    index = write_volume(tmp_path)
    label = Pds3Label(index)
    labels: list[Any] = []
    index_rows_ = GalileoSSI._index_rows

    def spy(filepath: Any, **kwargs: Any) -> Any:
        labels.append(kwargs['label'])
        return index_rows_(filepath, **kwargs)

    def parse(*args: Any) -> None:
        raise AssertionError('the index label was parsed again')

    monkeypatch.setattr(GalileoSSI, '_index_rows', spy)
    monkeypatch.setattr(oops.Host, '_parse_index_label', parse)
    GalileoSSI.from_index(index, label=label)
    assert labels == [label]


def test_from_index_rejects_an_unknown_keyword(tmp_path: pathlib.Path) -> None:
    with pytest.raises(TypeError, match="unexpected keyword argument 'supplementary'"):
        GalileoSSI.from_index(write_volume(tmp_path), supplementary=None)


@pytest.mark.parametrize('select', [0, slice(None), 'data'])
def test_from_index_refuses_select(tmp_path: pathlib.Path, select: Any) -> None:
    """The call itself raises, before any iteration."""

    with pytest.raises(oops.OopsValueError,
                       match='each row of a Galileo SSI index is one observation'):
        GalileoSSI.from_index(write_volume(tmp_path), select=select)

##########################################################################################
