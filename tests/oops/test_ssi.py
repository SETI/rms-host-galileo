##########################################################################################
# tests/oops/test_ssi.py
##########################################################################################
"""Tests for the Galileo SSI host, host_galileo.oops.ssi: the package, the detectors,
initialization, the clock-count start time, and the Snapshot that metadata give."""

import os
import pathlib
import subprocess
import sys
from typing import Any

import cspyce
import julian
import numpy as np
import oops
import pytest
from conftest import (EXPOSURE_MS, IMAGE_TIME, SCLK_COUNT, SSI_DISTORTION,
                      SSI_FOCAL_LENGTH, SSI_PIXEL_SIZE, SSI_PIXELS, SUMMED_PIXELS,
                      FakeSpice, FakeSpiceDB, fake_scs2e, image_array, label_dict)

import host_galileo
import host_galileo.oops
import host_galileo.oops.ssi as ssi
from host_galileo.oops._galileo import _Galileo
from host_galileo.oops.ssi      import GalileoSSI

OVERRIDE_NAMES = ('target', 'lightsource', 'path', 'frame', 'cadence', 'fov',
                  'calibrations', 'timeshift', 'navigation', 'tracker', 'parallel')


def make_snapshot(dict_: dict[str, Any], *, data: Any = None, full_fov: bool = False,
                  return_all_planets: bool = False,
                  **overrides: Any) -> oops.observation.Snapshot:
    """A Snapshot from `_make_snapshot`, with every override not given set to None."""

    GalileoSSI.initialize()
    obs: oops.observation.Snapshot = GalileoSSI._make_snapshot(
        dict_, GalileoSSI._interpret(dict_), data=data, full_fov=full_fov,
        return_all_planets=return_all_planets,
        overrides={name: overrides.get(name) for name in OVERRIDE_NAMES})
    return obs


def snapshot(**changes: Any) -> oops.observation.Snapshot:
    """A Snapshot made directly from the test label, without any file."""

    return make_snapshot(label_dict(**changes))


def tdb_from_image_time(image_time: str = IMAGE_TIME) -> float:
    """The time in seconds TDB of an IMAGE_TIME."""

    return float(julian.tdb_from_tai(julian.tai_from_iso(image_time)))


def tdb_from_count(count: str) -> float:
    """The time in seconds TDB of a label-form clock count, by the fake clock."""

    (rim, mod91) = count.split('.')
    return fake_scs2e(-77, f'{rim}:{mod91}:0:0')


def run_python(script: str) -> list[str]:
    """The lines that `script` prints when run by a new Python with the source on its
    path."""

    src = pathlib.Path(__file__).parents[2] / 'src'
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join([str(src), env.get('PYTHONPATH', '')])
    result = subprocess.run([sys.executable, '-c', script], env=env,
                            capture_output=True, text=True, check=True)
    return result.stdout.splitlines()

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
    assert run_python(_LAZY_IMPORT_SCRIPT) == ['False False', 'True True', 'True']


# Counts the kernels furnished by importing the host
_IMPORT_KERNELS_SCRIPT = """
import cspyce
furnished = []
cspyce.furnsh = furnished.append
import host_galileo.oops
print(len(furnished))
"""


def test_import_loads_no_kernels() -> None:
    assert run_python(_IMPORT_KERNELS_SCRIPT) == ['0']


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
    ({'INSTRUMENT_HOST_NAME': 'GALILEO ORBITER', 'INSTRUMENT_ID': 'SSI'}, True),
    ({'INSTRUMENT_HOST_NAME': 'GALILEO ORBITER', 'INSTRUMENT_ID': 'NIMS'}, False),
    ({'INSTRUMENT_HOST_NAME': 'CASSINI ORBITER', 'INSTRUMENT_ID': 'SSI'}, False),
], ids=['ssi', 'other-instrument', 'other-host'])
def test_detect_in_index(label: dict[str, Any], expected: bool) -> None:
    assert GalileoSSI._detect_in_index(label) is expected

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


def test_initialize_takes_keywords_only() -> None:
    with pytest.raises(TypeError, match='positional argument'):
        GalileoSSI.initialize([5])


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
# tstart_from_sclk()
##########################################################################################


@pytest.mark.parametrize(('frame_duration', 'delay'), [
    (2.333, 0.225), (8.667, 1.149), (15.166, 1.150), (30.333, 2.817), (60.667, 6.483),
])
def test_tstart_from_sclk_adds_the_shutter_delay(frame_duration: float,
                                                 delay: float) -> None:
    tstart = GalileoSSI.tstart_from_sclk(label_dict(FRAME_DURATION=frame_duration))
    assert tstart == pytest.approx(tdb_from_count(SCLK_COUNT) + delay, abs=1.e-6)


def test_tstart_from_sclk_reads_the_spice_form() -> None:
    row = label_dict(SPACECRAFT_CLOCK_START_COUNT='03496321:00:0:0')
    assert GalileoSSI.tstart_from_sclk(row) == pytest.approx(
        tdb_from_count(SCLK_COUNT) + 1.150, abs=1.e-6)


def test_tstart_from_sclk_opens_a_phase1_extended_exposure_a_frame_earlier() -> None:
    count = '03464058.90'
    tstart = GalileoSSI.tstart_from_sclk(label_dict(SPACECRAFT_CLOCK_START_COUNT=count,
                                                    FRAME_DURATION=30.333,
                                                    EXPOSURE_TYPE='EXTENDED'))
    assert tstart == pytest.approx(tdb_from_count(count) + 2.817 - 30.333, abs=1.e-6)


def test_tstart_from_sclk_opens_a_phase2_extended_exposure_at_the_normal_delay() -> None:
    count = '03464059.00'
    tstart = GalileoSSI.tstart_from_sclk(label_dict(SPACECRAFT_CLOCK_START_COUNT=count,
                                                    FRAME_DURATION=30.333,
                                                    EXPOSURE_TYPE='EXTENDED'))
    assert tstart == pytest.approx(tdb_from_count(count) + 2.817, abs=1.e-6)


def test_tstart_from_sclk_takes_a_missing_exposure_type_as_normal() -> None:
    count = '02492218.00'
    tstart = GalileoSSI.tstart_from_sclk(label_dict(SPACECRAFT_CLOCK_START_COUNT=count,
                                                    FRAME_DURATION=30.333,
                                                    EXPOSURE_TYPE=None))
    assert tstart == pytest.approx(tdb_from_count(count) + 2.817, abs=1.e-6)


def test_tstart_from_sclk_rejects_an_unknown_frame_duration() -> None:
    with pytest.raises(oops.host.HostFormatError,
                       match=r'unrecognized Galileo SSI FRAME_DURATION: 5\.0'):
        GalileoSSI.tstart_from_sclk(label_dict(FRAME_DURATION=5.))


def test_tstart_from_sclk_loads_only_the_clock_kernels(fake_spice: FakeSpice) -> None:
    GalileoSSI.tstart_from_sclk(label_dict())
    assert fake_spice.furnished == [_Galileo.SCLK_KERNEL]
    assert fake_spice.leap_seconds == 1
    assert not GalileoSSI._initialized

##########################################################################################
# The Snapshot that metadata give
##########################################################################################


def test_snapshot_time_is_centered_on_image_time() -> None:
    obs = snapshot()
    tmid = tdb_from_image_time()
    assert obs.time == pytest.approx((tmid - 0.13125, tmid + 0.13125))
    assert obs.texp == pytest.approx(EXPOSURE_MS / 1000.)
    assert obs.time_from_sclk is False


@pytest.mark.parametrize('image_time', [None, 'UNK'], ids=['absent', 'UNK'])
def test_snapshot_without_image_time_starts_from_the_clock(
        image_time: str | None) -> None:
    obs = snapshot(IMAGE_TIME=image_time)
    assert obs.time[0] == pytest.approx(tdb_from_count(SCLK_COUNT) + 1.150, abs=1.e-6)
    assert obs.time_from_sclk is True


def test_snapshot_keeps_a_zero_exposure() -> None:
    assert snapshot(EXPOSURE_DURATION=0.).texp == 0.


@pytest.mark.parametrize(('name', 'expected'), [
    ('host', 'Galileo'), ('instrument', 'SSI'), ('filter', 'VIOLET'),
    ('observation_id', 'G1GSGLOBAL02'), ('label_target', 'GANYMEDE'),
    ('spacecraft_clock_start_count', SCLK_COUNT), ('window', [1, 0, 400, 392]),
    ('full_fov', False), ('spice_frame_name', 'GLL_SCAN_PLATFORM'),
    ('spice_frame_id', -77001), ('ck_object_id', -77001), ('ck_sclk_id', -77),
    ('spice_frame_is_frozen', False),
])
def test_snapshot_subfield(name: str, expected: Any) -> None:
    assert getattr(snapshot(), name) == expected


def test_snapshot_dict_is_the_metadata() -> None:
    dict_ = label_dict()
    assert make_snapshot(dict_).dict is dict_


def test_snapshot_spice_to_frame_is_the_identity() -> None:
    assert snapshot().spice_to_frame == oops.Matrix3.IDENTITY


def test_snapshot_ck_ids_agree_with_spice() -> None:
    obs = snapshot()
    assert cspyce.ckmeta(obs.ck_object_id, 'SCLK') == obs.ck_sclk_id


@pytest.mark.parametrize(('item', 'name'), [
    ('FILTER_NAME', 'filter'), ('TARGET_NAME', 'label_target'),
    ('OBSERVATION_ID', 'observation_id'),
    ('SPACECRAFT_CLOCK_START_COUNT', 'spacecraft_clock_start_count'),
])
def test_snapshot_leaves_out_the_subfield_of_a_missing_item(item: str, name: str) -> None:
    assert name not in snapshot(**{item: None}).subfields


def test_snapshot_spice_kernels(fake_spicedb: FakeSpiceDB) -> None:
    assert snapshot().spice_kernels == fake_spicedb.used


def test_snapshot_geometry() -> None:
    obs = snapshot()
    assert obs.path.path_id == 'GLL'
    assert obs.frame.frame_id == 'GLL_SCAN_PLATFORM'
    assert obs.target == 'GANYMEDE'
    assert obs.lightsource == 'SUN'


def test_snapshot_uses_cks_for_the_default_frame(fake_spicedb: FakeSpiceDB) -> None:
    snapshot()
    assert fake_spicedb.named('used_basenames')[0][1]['types'] is None


def test_snapshot_uses_a_frame_override_as_given(fake_spice: FakeSpice,
                                                 fake_spicedb: FakeSpiceDB) -> None:
    # Unlike the fake scan-platform frame, this one is not an identity matrix
    frame = oops.frame.Cmatrix(oops.Matrix3([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),
                               oops.Frame.J2000, frame_id='CUSTOM')
    obs = make_snapshot(label_dict(), frame=frame)
    assert obs.frame is frame
    assert 'CK' not in fake_spicedb.named('used_basenames')[0][1]['types']
    assert fake_spice.type1_frames == []


def test_snapshot_return_all_planets(fake_spicedb: FakeSpiceDB) -> None:
    make_snapshot(label_dict(), return_all_planets=True)
    assert 1 in fake_spicedb.named('used_basenames')[0][1]['bodies']


def test_snapshot_cutout_window_trims_the_fov() -> None:
    obs = snapshot(TELEMETRY_FORMAT_ID='IM8')
    assert type(obs.fov) is oops.fov.SliceFOV
    assert tuple(np.asarray(obs.fov.uv_shape.vals)) == (392, 399)
    assert tuple(np.asarray(obs.fov.uv_los.vals)) == pytest.approx((400., 399.))


@pytest.mark.parametrize('window', [None, [-1, -1, -1, -1]], ids=['absent', 'no-window'])
def test_snapshot_without_cutout_window(window: list[int] | None) -> None:
    obs = snapshot(CUT_OUT_WINDOW=window)
    assert obs.fov is GalileoSSI._FOVS['HIS']
    assert obs.window == [0, 0, SUMMED_PIXELS, SUMMED_PIXELS]


def test_snapshot_full_fov_keeps_the_full_fov_and_gives_the_window() -> None:
    obs = make_snapshot(label_dict(), full_fov=True)
    assert obs.fov is GalileoSSI._FOVS['HIS']
    assert obs.window == [1, 0, 400, 392]
    assert obs.full_fov is True


def test_snapshot_trims_the_data() -> None:
    obs = make_snapshot(label_dict(), data=image_array())
    assert np.array_equal(obs.data, image_array()[1:400, 0:392])


def test_snapshot_full_fov_gives_the_summed_image() -> None:
    obs = make_snapshot(label_dict(), data=image_array(), full_fov=True)
    assert np.array_equal(obs.data, image_array()[:SUMMED_PIXELS, :SUMMED_PIXELS])


@pytest.mark.parametrize(('label_target', 'target'), [
    ('BLACK_SKY', 'NONE'), ('BLACK SKY', 'NONE'), ('NON_SCIENCE', 'NONE'),
    ('NON SCIENCE', 'NONE'), ('STAR', 'NONE'), ('J RINGS', 'JUPITER_RING_PLANE'),
    ('IO', 'IO'),
])
def test_snapshot_target(label_target: str, target: str) -> None:
    obs = snapshot(TARGET_NAME=label_target)
    assert obs.target == target
    assert obs.label_target == label_target
    assert obs.lightsource == 'SUN'


def test_snapshot_without_target() -> None:
    assert snapshot(TARGET_NAME=None).target is None


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
    obs = make_snapshot(label_dict(), tracker='midtime')
    assert type(obs.frame) is oops.frame.TrackerFrame


def test_snapshot_navigation() -> None:
    obs = make_snapshot(label_dict(), navigation=(0., 0.))
    assert type(obs.frame) is oops.frame.Navigation


def test_snapshot_timeshift() -> None:
    obs = make_snapshot(label_dict(), timeshift=(5., 'T'))
    assert type(obs.cadence) is oops.cadence.TimeShift

##########################################################################################
