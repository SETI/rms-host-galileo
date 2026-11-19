##########################################################################################
# tests/oops/test_ssi.py
##########################################################################################
"""Tests for the Galileo SSI host, host_galileo.oops.ssi."""

from collections.abc import Iterator
import json
import os
import pathlib
import subprocess
import sys
from typing import Any, cast

import julian
import numpy as np
import oops
from filecache import FCPath
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


_VICAR_LABEL: dict[str, int | float | str | list[int]] = {
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


def make_snapshot(dict_: dict[str, Any], *, filepath: str = 'X.LBL', data: Any = None,
                  full_fov: bool = False, return_all_planets: bool = False,
                  default_mode: str | None = None,
                  **overrides: Any) -> oops.observation.Snapshot:
    """A Snapshot from `_make_snapshot`, with every override not given set to None."""

    GalileoSSI.initialize()
    obs: oops.observation.Snapshot = GalileoSSI._make_snapshot(
        dict_, filepath=filepath, data=data, full_fov=full_fov,
        return_all_planets=return_all_planets, default_mode=default_mode,
        overrides={name: overrides.get(name) for name in OVERRIDE_NAMES})
    return obs


def snapshot(**changes: Any) -> oops.observation.Snapshot:
    """A Snapshot made from the test label as an index row, without any file."""

    return make_snapshot(_ssi_labels.row_items(label_dict(**changes)),
                         filepath='GO_0017/X.LBL')


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
    ({'IMAGE_INDEX_TABLE': {'ROWS': 1}}, False),
], ids=['ssi', 'other-instrument', 'other-host', 'unnamed'])
def test_detect_in_index(label: dict[str, Any], expected: bool) -> None:
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
    assert obs.filepath == 'GO_0017/X.LBL'
    assert obs.basename == 'X.LBL'


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
    obs = snapshot(CUT_OUT_WINDOW=np.array([-1, -1, -1, -1]))
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


def test_snapshot_index_null_is_no_value() -> None:
    obs = snapshot(IMAGE_TIME='UNK')
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
                             r'GO_0017/X\.LBL'):
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


def test_from_file_sets_local_paths(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path), astrometry=True)
    image = tmp_path / 'C0349632100R.IMG'
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


def test_local_label_astrometry_without_the_image_names_it(
        tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path)
    (tmp_path / 'C0349632100R.IMG').unlink()
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert obs.abspath == (tmp_path / 'C0349632100R.IMG').resolve()


# Treated as remote, the image file named by the label is a file not in the cache unless
# the full read fetched it or an earlier read left it there
@pytest.mark.parametrize(('astrometry', 'cached', 'has_abspath'),
                         [(True, False, False), (True, True, True), (False, True, True)],
                         ids=['astrometry', 'astrometry-cached', 'full'])
def test_remote_label_abspath(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch,
                              astrometry: bool, cached: bool, has_abspath: bool) -> None:
    path = write_pds3_image(tmp_path)
    if not cached:
        (tmp_path / 'C0349632100R.IMG').unlink()
    monkeypatch.setattr(FCPath, 'is_local', lambda self: False)
    obs = GalileoSSI.from_file(path, astrometry=astrometry)
    assert hasattr(obs, 'abspath') is has_abspath


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
    assert (info.value.spice_code, str(info.value.filepath)) == ('SPICE(CKINSUFFDATA)',
                                                                 str(path))
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
                       match=r'no PDS3 parser can read .*C0349632100R\.LBL'):
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

_SUPPLEMENT_COLUMNS = [
    ('VOLUME_ID', 'CHARACTER', 7),
    ('FILE_SPECIFICATION_NAME', 'CHARACTER', 30),
    ('TELEMETRY_FORMAT_ID', 'CHARACTER', 3),
    ('CUT_OUT_WINDOW', 'ASCII_INTEGER', 4),
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
                columns: list[tuple[str, str, int]], rows: list[tuple[Any, ...]],
                **keywords: Any) -> pathlib.Path:
    """Write an ASCII PDS3 table and its label; return the label's path.

    An ASCII_INTEGER column holds a list of four values, each in a field of `width`
    bytes; an ASCII_REAL column holds one value.
    """

    def fields(value: Any, kind: str, width: int) -> list[str]:
        if kind == 'CHARACTER':
            return ['"' + str(value).ljust(width) + '"']
        if kind == 'ASCII_INTEGER':
            return [f'{item:{width}d}' for item in value]
        return [f'{value:{width}.3f}']

    lines = [','.join(item for (value, (_, kind, width)) in zip(row, columns, strict=True)
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
    for (column, kind, width) in columns:
        text += ['  OBJECT = COLUMN', f'    NAME = {column}', f'    DATA_TYPE = {kind}']
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


def write_index(directory: pathlib.Path) -> pathlib.Path:
    """Write a three-row SSI index and its label; return the label's path."""

    return write_table(directory, 'GO_0017_index', 'IMAGE_INDEX_TABLE', _INDEX_COLUMNS,
                       _INDEX_ROWS, **_INDEX_KEYWORDS)


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

    rows = GalileoSSI._read_index_rows(write_index(tmp_path))
    del rows[1]['EXPOSURE_DURATION']
    snapshots = GalileoSSI.from_index(tmp_path / 'missing.lbl', row_dicts=rows)
    assert next(snapshots).basename == 'C0346405900R.IMG'
    with pytest.raises(oops.host.HostFormatError,
                       match=r'exposure unknown: no EXPOSURE_DURATION for '
                             r'.*C0349632100R\.IMG'):
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
    obs = list(GalileoSSI.from_index(write_volume(tmp_path)))[1]
    image = tmp_path / 'GO_0017/G1/GANYMEDE/C0349632100R.IMG'
    assert (obs.basename, obs.image_url, obs.abspath) == ('C0349632100R.IMG',
                                                          image.as_posix(),
                                                          image.resolve())


# Treated as remote, the image file is not in the cache, since an index read fetches none
def test_from_index_remote_image_has_no_abspath(tmp_path: pathlib.Path,
                                                monkeypatch: pytest.MonkeyPatch) -> None:
    index = write_volume(tmp_path)
    monkeypatch.setattr(FCPath, 'is_local', lambda self: False)
    assert not hasattr(list(GalileoSSI.from_index(index))[1], 'abspath')


@pytest.mark.parametrize(('index', 'volumes'), [
    ('/h/volumes/GO_0xxx/GO_0017/INDEX/IMGINDEX.LBL', '/h/volumes/GO_0xxx'),
    ('/h/metadata/GO_0xxx/GO_0017/GO_0017_index.lbl', '/h/volumes/GO_0xxx'),
    ('https://h.org/metadata/GO_0xxx/GO_0017/GO_0017_index.lbl',
     'https://h.org/volumes/GO_0xxx'),
    ('/h/GO_0017_index.lbl', '/h'),
], ids=['volume-index', 'metadata', 'metadata-url', 'other'])
def test_volumes_directory(index: str, volumes: str) -> None:
    assert str(GalileoSSI._volumes_directory(index)) == volumes


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


def test_from_index_rejects_a_mismatched_supplement(tmp_path: pathlib.Path) -> None:
    rows = [_SUPPLEMENT_ROWS[1], _SUPPLEMENT_ROWS[0], _SUPPLEMENT_ROWS[2]]
    with pytest.raises(oops.host.HostIndexError,
                       match=r'FILE_SPECIFICATION_NAME value mismatch.*\n'
                             r'index file mismatch at record \[0\]'):
        GalileoSSI.from_index(write_index(tmp_path),
                              supplement=write_supplement(tmp_path, rows))


def test_from_index_rejects_a_short_supplement(tmp_path: pathlib.Path) -> None:
    with pytest.raises(oops.host.HostIndexError,
                       match='index files have different record counts'):
        GalileoSSI.from_index(write_index(tmp_path),
                              supplement=write_supplement(tmp_path, _SUPPLEMENT_ROWS[:2]))


def test_from_index_filter(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_volume(tmp_path),
                                      filter={'FILTER_NAME': 'VIOLET'})
    assert [obs.filter for obs in snapshots] == ['VIOLET']


def test_from_index_uses_given_rows(tmp_path: pathlib.Path) -> None:
    rows = GalileoSSI._read_index_rows(write_index(tmp_path))
    snapshots = GalileoSSI.from_index(tmp_path / 'missing.lbl', row_dicts=rows[:1])
    assert [obs.basename for obs in snapshots] == ['C0346405900R.IMG']


def test_from_index_names_the_image_in_the_label_case(tmp_path: pathlib.Path) -> None:
    rows = GalileoSSI._read_index_rows(write_index(tmp_path))
    rows[0]['FILE_SPECIFICATION_NAME'] = 'j0/opnav/c0346405900r.lbl'
    [obs] = GalileoSSI.from_index(tmp_path / 'missing.lbl', row_dicts=rows[:1])
    assert obs.basename == 'c0346405900r.img'


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
    snapshots = oops.Host.from_index(write_index(tmp_path))
    assert [obs.instrument for obs in snapshots] == ['SSI', 'SSI', 'SSI']


# The supplement is required, so None is refused rather than taken as "search for it"
def test_from_index_refuses_no_supplement(tmp_path: pathlib.Path) -> None:
    with pytest.raises(oops.host.HostIndexError,
                       match=r'Galileo SSI index needs its supplemental index: '
                             r'.*GO_0017_index\.lbl'):
        GalileoSSI.from_index(write_volume(tmp_path), supplement=None)


def test_from_index_rows_need_no_supplement(tmp_path: pathlib.Path) -> None:
    rows = GalileoSSI._read_index_rows(write_index(tmp_path))
    assert len(list(GalileoSSI.from_index(tmp_path / 'x.lbl', row_dicts=rows,
                                          supplement=None))) == 3


def test_from_index_uses_a_given_label(tmp_path: pathlib.Path,
                                       monkeypatch: pytest.MonkeyPatch) -> None:
    index = write_volume(tmp_path)
    label = Pds3Label(index)
    labels: list[Any] = []
    read_index_rows = GalileoSSI._read_index_rows

    def spy(filepath: Any, **kwargs: Any) -> Any:
        labels.append(kwargs['label'])
        return read_index_rows(filepath, **kwargs)

    monkeypatch.setattr(GalileoSSI, '_read_index_rows', spy)
    GalileoSSI.from_index(index, label=label)
    assert labels == [label]


def test_from_index_rejects_an_unknown_keyword(tmp_path: pathlib.Path) -> None:
    with pytest.raises(oops.OopsTypeError,
                       match=r"unexpected keyword argument for GalileoSSI\.from_index: "
                             r"'supplementary'"):
        GalileoSSI.from_index(write_volume(tmp_path), supplementary=None)


@pytest.mark.parametrize('select', [0, slice(None), 'data'])
def test_from_index_refuses_select(tmp_path: pathlib.Path, select: Any) -> None:
    """The call itself raises, before any iteration."""

    with pytest.raises(oops.OopsValueError,
                       match='each row of a Galileo SSI index is one observation'):
        GalileoSSI.from_index(write_volume(tmp_path), select=select)

##########################################################################################
