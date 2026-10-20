##########################################################################################
# tests/oops/test_ssi.py
##########################################################################################
"""Tests for the Galileo SSI host, host_galileo.oops.ssi."""

import datetime as dt
import os
import pathlib
import subprocess
import sys
from typing import Any, cast

import julian
import numpy as np
import oops
import pytest
import vicar
from conftest import (SSI_DISTORTION, SSI_FOCAL_LENGTH, SSI_PIXEL_SIZE, SSI_PIXELS,
                      FakeSpice, FakeSpiceDB, fake_scs2e, pds3_lines, write_text)

import host_galileo
import host_galileo.oops
import host_galileo.oops.ssi as ssi
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


def write_pds3_label(directory: pathlib.Path, stem: str, **changes: Any) -> pathlib.Path:
    """Write a detached PDS3 label for the 800x800 image `stem`.IMG; return its path."""

    lines = ['PDS_VERSION_ID = PDS3', 'RECORD_TYPE = FIXED_LENGTH',
             f'RECORD_BYTES = {SSI_PIXELS}', f'FILE_RECORDS = {SSI_PIXELS}',
             f'^IMAGE = "{stem}.IMG"']
    lines += pds3_lines(label_dict(**changes))
    lines += ['OBJECT = IMAGE', f'  LINES = {SSI_PIXELS}',
              f'  LINE_SAMPLES = {SSI_PIXELS}', '  SAMPLE_BITS = 8',
              '  SAMPLE_TYPE = UNSIGNED_INTEGER', 'END_OBJECT = IMAGE']
    return write_text(directory / f'{stem}.LBL', lines)


def write_pds3_image(directory: pathlib.Path, **changes: Any) -> pathlib.Path:
    """Write an image as raw bytes with its detached PDS3 label; return the label path."""

    (directory / 'C0349632100R.IMG').write_bytes(image_array().tobytes())
    return write_pds3_label(directory, 'C0349632100R', **changes)


_VICAR_LABEL: dict[str, int | float | str | list[int]] = {
    'MISSION'       : 'GALILEO',
    'SENSOR'        : 'SSI',
    'RIM'           : 3496321,
    'MOD91'         : 0,
    'PA'            : 'G1GSGLOBAL02',
    'TARGET'        : 'GANYMEDE',
    'FILTER'        : 3,
    'EXP'           : EXPOSURE_MS,
    'TLMFMT'        : 'HIS',
    'CUT_OUT_WINDOW': [2, 1, 399, 392],
    'SCETYEAR'      : 1996,
    'SCETDAY'       : 178,
    'SCETHOUR'      : 8,
    'SCETMIN'       : 46,
    'SCETSEC'       : 4,
    'SCETMSEC'      : 890,
}


def write_vicar_image(directory: pathlib.Path, **extra: int) -> pathlib.Path:
    """Write an image as a VICAR file with no PDS3 label; return the image path.

    The header holds the VICAR keywords of C0349632100R that the host reads, followed by
    any `extra` keywords.
    """

    vic = vicar.VicarImage.from_array(image_array())
    for (key, value) in {**_VICAR_LABEL, **extra}.items():
        vic[key] = value
    path = directory / 'C0349632100R.IMG'
    vic.write_file(path)
    return path


OVERRIDE_NAMES = ('target', 'lightsource', 'path', 'frame', 'cadence', 'fov',
                  'calibrations', 'timeshift', 'navigation', 'tracker', 'parallel')


def make_snapshot(dict_: dict[str, Any], *, filepath: str = 'X.LBL', data: Any = None,
                  full_fov: bool = False, return_all_planets: bool = False,
                  **overrides: Any) -> oops.observation.Snapshot:
    """A Snapshot from `_make_snapshot`, with every override not given set to None."""

    GalileoSSI.initialize()
    obs: oops.observation.Snapshot = GalileoSSI._make_snapshot(
        dict_, filepath=filepath, data=data, full_fov=full_fov,
        return_all_planets=return_all_planets,
        overrides={name: overrides.get(name) for name in OVERRIDE_NAMES})
    return obs


def snapshot(**changes: Any) -> oops.observation.Snapshot:
    """A Snapshot made directly from the test label, without any file."""

    return make_snapshot(label_dict(**changes), filepath='GO_0017/X.LBL')


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
    ({'SPACECRAFT_NAME': 'GALILEO ORBITER', 'INSTRUMENT_NAME': 'SOLID STATE IMAGING SYSTEM'},
     True),
    ({'SPACECRAFT_NAME': 'GALILEO ORBITER', 'INSTRUMENT_NAME': 'SOLID_STATE_IMAGING'}, True),
    ({'SPACECRAFT_NAME': 'GALILEO ORBITER', 'INSTRUMENT_NAME': 'NEAR INFRARED MAPPING'},
     False),
    ({'SPACECRAFT_NAME': 'CASSINI ORBITER', 'INSTRUMENT_NAME': 'SOLID STATE IMAGING'}, False),
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
    assert GalileoSSI._INSTRUMENT_KERNEL is fake_spicedb.kernel_dict
    assert fake_spice.furnished[0] == 'Galileo/SCLK/mk00062a.tsc'
    assert GalileoSSI._initialized


def test_initialize_rejects_an_unexpected_image_shape(fake_spicedb: FakeSpiceDB) -> None:
    fake_spicedb.kernel_dict['INS'][-77036]['MAX_LINE'] = 400
    with pytest.raises(oops.host.HostError,
                       match=r'unsupported SSI image shape \(800, 400\) in the SSI kernel'):
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
    assert list(fov.coefft_uv_from_xy) == pytest.approx([1., 0., SSI_DISTORTION / scale**2])


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


def test_snapshot_missing_telemetry_format() -> None:
    obs = snapshot(TELEMETRY_FORMAT_ID=None, CUT_OUT_WINDOW=None)
    assert obs.fov is GalileoSSI._FOVS['NONE']


def test_snapshot_trims_the_data() -> None:
    GalileoSSI.initialize()
    obs = make_snapshot(label_dict(), data=image_array())
    assert obs.data.shape == (399, 392)
    assert obs.data[0, 0] == image_array()[1, 0]


@pytest.mark.parametrize(('label_target', 'target'), [
    ('GANYMEDE', 'GANYMEDE'),
    ('BLACK_SKY', 'NONE'),
    ('NON_SCIENCE', 'NONE'),
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
    assert obs.time[0] == pytest.approx(fake_scs2e(-77, '3496321:0:0:0') + delay)


def test_extended_exposure_starts_a_frame_earlier() -> None:
    normal = snapshot(IMAGE_TIME='UNK')
    extended = snapshot(IMAGE_TIME='UNK', EXPOSURE_TYPE='EXTENDED')
    assert normal.time[0] - extended.time[0] == pytest.approx(15.167)


def test_missing_exposure_type_is_normal() -> None:
    normal = snapshot(IMAGE_TIME='UNK')
    missing = snapshot(IMAGE_TIME='UNK', EXPOSURE_TYPE=None)
    assert missing.time == normal.time


def test_unknown_frame_duration_is_rejected() -> None:
    with pytest.raises(oops.host.HostFormatError,
                       match=r'unrecognized Galileo SSI FRAME_DURATION: 20\.0'):
        snapshot(IMAGE_TIME='UNK', FRAME_DURATION=20.)

##########################################################################################
# from_file()
##########################################################################################


def test_from_file_pds3(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path))
    assert type(obs) is oops.observation.Snapshot
    assert obs.data.shape == (399, 392)
    assert np.array_equal(obs.data, image_array()[1:400, 0:392])
    assert obs.basename == 'C0349632100R.LBL'
    assert obs.observation_id == 'G1GSGLOBAL02'


def test_from_file_reads_the_window_in_label_order(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path, CUT_OUT_WINDOW=[1, 2, 300, 200]))
    assert obs.data.shape == (300, 200)
    assert obs.data[0, 0] == image_array()[0, 1]


def test_from_file_full_fov(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_pds3_image(tmp_path), full_fov=True)
    assert obs.data.shape == (SSI_PIXELS, SSI_PIXELS)
    assert obs.fov is GalileoSSI._FOVS['HIS']


def test_from_file_sets_local_paths(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path)
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert obs.abspath == path.resolve()
    assert obs.image_url == path.as_posix()


def test_from_file_pds3_astrometry(tmp_path: pathlib.Path) -> None:
    path = write_pds3_image(tmp_path)
    (tmp_path / 'C0349632100R.IMG').unlink()
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert not hasattr(obs, 'data')
    assert obs.filter == 'VIOLET'


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
    assert obs.dict['EARTH_RECEIVED_TIME'] == dt.datetime(1996, 6, 26, 9, 30, 15, 250000)


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


def test_from_file_accepts_a_parsed_label(tmp_path: pathlib.Path) -> None:
    from pdsparser import Pds3Label
    label = Pds3Label(write_pds3_image(tmp_path))
    assert GalileoSSI.from_file(label).data.shape == (399, 392)


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
    with pytest.raises(oops.host.UnknownHostError, match='unrecognized host'):
        GalileoSSI.from_file(path, astrometry=True)

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
                     '    ITEMS = 4', f'    ITEM_BYTES = {width}', f'    ITEM_OFFSET = {width + 1}']
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

    return write_table(directory, 'GO_0017_supplemental_index', 'SUPPLEMENTAL_INDEX_TABLE',
                       _SUPPLEMENT_COLUMNS, rows, **_INDEX_KEYWORDS)


def test_from_index_keeps_zero_exposures(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_index(tmp_path))
    assert [obs.basename for obs in snapshots] == ['C0346405900R.LBL', 'C0349632100R.LBL',
                                                   'C0349632200R.LBL']
    assert snapshots[2].texp == 0.


def test_from_index_rows(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_index(tmp_path))
    assert [obs.filepath for obs in snapshots] == ['GO_0017/J0/OPNAV/C0346405900R.LBL',
                                                   'GO_0017/G1/GANYMEDE/C0349632100R.LBL',
                                                   'GO_0017/G1/GANYMEDE/C0349632200R.LBL']
    assert [obs.filter for obs in snapshots] == ['CLEAR', 'VIOLET', 'CLEAR']
    assert [obs.observation_id for obs in snapshots] == ['J0GSOPNAV115', 'G1GSGLOBAL02',
                                                         'G1GSGLOBAL02']
    assert snapshots[1].time[0] == pytest.approx(tdb_from_image_time() - 0.13125)


def test_from_index_without_supplement_uses_the_full_fov(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_index(tmp_path))
    assert all(obs.fov is GalileoSSI._FOVS['NONE'] for obs in snapshots)


def test_from_index_with_supplement(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_index(tmp_path),
                                      supplement=write_supplement(tmp_path))
    assert snapshots[0].fov is GalileoSSI._FOVS['IM8']
    assert tuple(snapshots[1].fov.uv_shape.vals) == (392, 399)


def test_from_index_keeps_the_index_values(tmp_path: pathlib.Path) -> None:
    snapshots = GalileoSSI.from_index(write_index(tmp_path),
                                      supplement=write_supplement(tmp_path))
    assert snapshots[1].dict['SPACECRAFT_CLOCK_START_COUNT'] == SCLK_COUNT


def test_from_index_finds_the_supplement(tmp_path: pathlib.Path) -> None:
    write_supplement(tmp_path)
    snapshots = GalileoSSI.from_index(write_index(tmp_path))
    assert snapshots[0].fov is GalileoSSI._FOVS['IM8']


def test_from_index_requires_a_named_supplement(tmp_path: pathlib.Path) -> None:
    with pytest.raises(FileNotFoundError, match=r'missing_supplemental_index\.lbl'):
        GalileoSSI.from_index(write_index(tmp_path),
                              supplement=tmp_path / 'missing_supplemental_index.lbl')


def test_from_index_rejects_another_host(tmp_path: pathlib.Path) -> None:
    path = write_table(tmp_path, 'GO_0017_index', 'IMAGE_INDEX_TABLE', _INDEX_COLUMNS,
                       _INDEX_ROWS, MISSION_NAME='CASSINI', INSTRUMENT_HOST_NAME='CASSINI',
                       INSTRUMENT_ID='ISS')
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
    snapshots = GalileoSSI.from_index(write_index(tmp_path),
                                      filter={'FILTER_NAME': 'VIOLET'})
    assert [obs.filter for obs in snapshots] == ['VIOLET']


def test_from_index_uses_given_rows(tmp_path: pathlib.Path) -> None:
    rows = GalileoSSI._read_index_rows(write_index(tmp_path))
    snapshots = GalileoSSI.from_index(tmp_path / 'missing.lbl', row_dicts=rows[:1])
    assert [obs.basename for obs in snapshots] == ['C0346405900R.LBL']


def test_from_index_fov_override(tmp_path: pathlib.Path) -> None:
    fov = oops.fov.FlatFOV((1.e-5, 1.e-5), (10, 10))
    snapshots = GalileoSSI.from_index(write_index(tmp_path), fov=fov)
    assert all(obs.fov is fov for obs in snapshots)


def test_from_index_parallel(tmp_path: pathlib.Path) -> None:
    parallel = snapshot()
    snapshots = GalileoSSI.from_index(write_index(tmp_path), parallel=parallel)
    assert all(type(obs.frame) is oops.frame.Cmatrix for obs in snapshots)
    assert all(obs.frame.reference is parallel.frame for obs in snapshots)


def test_host_from_index_dispatches_to_ssi(tmp_path: pathlib.Path) -> None:
    snapshots = oops.Host.from_index(write_index(tmp_path))
    assert [obs.instrument for obs in snapshots] == ['SSI', 'SSI', 'SSI']

##########################################################################################
