##########################################################################################
# tests/oops/conftest.py
##########################################################################################
"""Shared fixtures that let the Galileo host run without SPICE kernels or spicedb.

Every spicedb function the package calls is replaced by :class:`FakeSpiceDB`, and the
oops and cspyce calls that would read SPICE kernels are replaced by kernel-free
equivalents. Everything else, including oops itself, runs for real.
"""

import pathlib
import re
from dataclasses import dataclass, field
from typing import Any

import cspyce
import julian
import oops
import pytest
import spicedb

import host_galileo.oops._galileo as galileo_module
from host_galileo.oops.ssi import GalileoSSI

# Values from the SSI instrument kernel, gll36001.ti
SSI_DISTORTION = 6.58e-9
SSI_FOCAL_LENGTH = 1501.039         # mm
SSI_PIXEL_SIZE = 15.24e-3           # mm
SSI_PIXELS = 800

# The fake spacecraft clock: one RIM count is one 60.667 s frame of the mod-91 counter,
# starting at SCLK_EPOCH
SCLK_EPOCH = julian.tdb_from_iso('1989-10-18T00:00:00')
RIM_SECONDS = 60.667

_ISO_TIME = re.compile(r'\d{4}-\d{2}-\d{2}T[\d:.]+Z?')


def pds3_value(value: Any) -> str:
    """A Python value formatted for a PDS3 label; ISO times and numbers stay unquoted.

    A tuple becomes a PDS3 sequence and a list a PDS3 set, which is how labels write
    CUT_OUT_WINDOW.
    """

    if isinstance(value, tuple):
        return '(' + ', '.join(pds3_value(item) for item in value) + ')'
    if isinstance(value, list):
        return '{' + ', '.join(pds3_value(item) for item in value) + '}'
    if isinstance(value, str) and not _ISO_TIME.fullmatch(value):
        return f'"{value}"'
    return str(value)


def pds3_lines(keywords: dict[str, Any], indent: str = '') -> list[str]:
    """PDS3 label lines assigning each keyword its formatted value."""

    return [f'{indent}{key} = {pds3_value(value)}' for (key, value) in keywords.items()]


def fake_scs2e(sc: int, count: str) -> float:
    """The fake spacecraft clock conversion; fails unless `sc` is the Galileo ID."""

    assert sc == -77, f'unexpected spacecraft ID {sc}'
    (rim, mod91, mod10, mod8) = (int(field) for field in count.split(':'))
    return float(SCLK_EPOCH + RIM_SECONDS * (rim + (mod91 + (mod10 + mod8/8) / 10) / 91))


def ik_dict() -> dict[str, Any]:
    """The part of the SSI instrument kernel that the host reads as a dictionary."""

    return {'INS': {-77036: {'MAX_SAMPLE': SSI_PIXELS, 'MAX_LINE': SSI_PIXELS}}}


_GDPOOL = {
    'INS-77036_DISTORTION_COEFF': [SSI_DISTORTION],
    'INS-77036_FOCAL_LENGTH'    : [SSI_FOCAL_LENGTH],
    'INS-77036_PIXEL_SIZE'      : [SSI_PIXEL_SIZE],
    'INS-77036_FOV_CENTER'      : [400., 400.],
}


@dataclass
class FakeKernel:
    """A stand-in for a spicedb KernelInfo object."""

    filespec: str
    start_time: str = '1995-01-01T00:00:00.000'
    stop_time: str = '1995-01-01T00:00:00.000'


@dataclass
class FakeSpiceDB:
    """Records every spicedb call and returns configurable results.

    Attributes:
        calls (list[tuple[str, tuple, dict]]): Each call as (name, args, kwargs).
        inst_kernels (list[FakeKernel]): Returned by `select_inst`.
        kernel_dict (dict): Returned by `as_dict`.
        used (list[str]): Returned by `used_basenames`.
    """

    calls: list[tuple[str, tuple[Any, ...], dict[str, Any]]] = field(default_factory=list)
    inst_kernels: list[FakeKernel] = field(
        default_factory=lambda: [FakeKernel('gll36001.ti')])
    kernel_dict: dict[str, Any] = field(default_factory=ik_dict)
    used: list[str] = field(default_factory=lambda: ['gll36001.ti', 'mk00062a.tsc'])

    def _record(self, func: str, /, *args: Any, **kwargs: Any) -> None:
        self.calls.append((func, args, kwargs))

    def named(self, name: str) -> list[tuple[tuple[Any, ...], dict[str, Any]]]:
        """The (args, kwargs) of every recorded call to the named function."""

        return [(args, kwargs) for (n, args, kwargs) in self.calls if n == name]

    def open_db(self, *args: Any, **kwargs: Any) -> None:
        self._record('open_db', *args, **kwargs)

    def close_db(self, *args: Any, **kwargs: Any) -> None:
        self._record('close_db', *args, **kwargs)

    def select_inst(self, *args: Any, **kwargs: Any) -> list[FakeKernel]:
        self._record('select_inst', *args, **kwargs)
        return self.inst_kernels

    def furnish_inst(self, *args: Any, **kwargs: Any) -> list[str]:
        self._record('furnish_inst', *args, **kwargs)
        return []

    def furnish_kernels(self, *args: Any, **kwargs: Any) -> list[str]:
        self._record('furnish_kernels', *args, **kwargs)
        return []

    def as_dict(self, *args: Any, **kwargs: Any) -> dict[str, Any]:
        self._record('as_dict', *args, **kwargs)
        return self.kernel_dict

    def as_names(self, kernels: list[FakeKernel]) -> list[str]:
        self._record('as_names', kernels)
        return [kernel.filespec for kernel in kernels]

    def used_basenames(self, *args: Any, **kwargs: Any) -> list[str]:
        self._record('used_basenames', *args, **kwargs)
        return self.used


_SPICEDB_FUNCTIONS = ('open_db', 'close_db', 'select_inst', 'furnish_inst',
                      'furnish_kernels', 'as_dict', 'as_names', 'used_basenames')


class FakeFileCachePrefix:
    """A stand-in for the SPICE FCPath prefix; `retrieve` returns paths unchanged."""

    def __init__(self, retrieved: list[Any]) -> None:
        self._retrieved = retrieved

    def retrieve(self, paths: str | list[str]) -> str | list[str]:
        self._retrieved.append(paths)
        return paths


@dataclass
class FakeSpice:
    """Records the oops and cspyce calls that would otherwise read SPICE kernels.

    Attributes:
        leap_seconds (int): The number of calls to `oops.spice.load_leap_seconds`.
        solar_system (list[tuple[tuple, dict]]): The (args, kwargs) of each call to
            `Body.define_solar_system`.
        spice_paths (list[tuple[str, str]]): The (path ID, origin) of each SpicePath.
        type1_frames (list[tuple[str, int]]): The (frame ID, tick tolerance) of each
            SpiceType1Frame.
        retrieved (list): The argument of each call to the SPICE prefix's `retrieve`.
        furnished (list[str]): The path given to each call to `cspyce.furnsh`.
    """

    leap_seconds: int = 0
    solar_system: list[tuple[tuple[Any, ...], dict[str, Any]]] = field(
        default_factory=list)
    spice_paths: list[tuple[str, str]] = field(default_factory=list)
    type1_frames: list[tuple[str, int]] = field(default_factory=list)
    retrieved: list[Any] = field(default_factory=list)
    furnished: list[str] = field(default_factory=list)


@pytest.fixture(autouse=True)
def fake_spicedb(monkeypatch: pytest.MonkeyPatch) -> FakeSpiceDB:
    """Replace every spicedb function the package calls with a recording fake."""

    fake = FakeSpiceDB()
    for name in _SPICEDB_FUNCTIONS:
        monkeypatch.setattr(spicedb, name, getattr(fake, name))
    return fake


@pytest.fixture(autouse=True)
def fake_spice(monkeypatch: pytest.MonkeyPatch) -> FakeSpice:
    """Keep SPICE kernels out of every test and isolate the oops registries.

    Kernel-reading calls are replaced: time conversion uses julian, the clock conversion
    uses :func:`fake_scs2e`, the instrument kernel pool holds the SSI values, a SpicePath
    becomes a FixedPath at the origin, and a SpiceType1Frame becomes an identity Cmatrix
    on J2000. The SPICE file cache returns paths without reading them, and
    `cspyce.furnsh` only records the path. The Path and Frame registries are copied so
    that frames registered by one test never reach another.
    """

    fake = FakeSpice()

    def load_leap_seconds() -> None:
        fake.leap_seconds += 1

    def define_solar_system(*args: Any, **kwargs: Any) -> None:
        fake.solar_system.append((args, kwargs))

    def spice_path(spice_id: str, origin: str) -> oops.Path:
        fake.spice_paths.append((spice_id, origin))
        return oops.path.FixedPath((0., 0., 0.), oops.Path.SSB, oops.Frame.J2000,
                                   path_id=spice_id)

    def spice_type1_frame(frame_id: str, tick_tolerance: int) -> oops.Frame:
        fake.type1_frames.append((frame_id, tick_tolerance))
        return oops.frame.Cmatrix(oops.Matrix3.IDENTITY, oops.Frame.J2000, frame_id=frame_id)

    def gdpool(name: str, start: int) -> list[float]:
        return _GDPOOL[name][start:]

    def furnsh(path: str) -> None:
        fake.furnished.append(path)

    monkeypatch.setattr(cspyce, 'str2et', julian.tdb_from_iso)
    monkeypatch.setattr(cspyce, 'scs2e', fake_scs2e)
    monkeypatch.setattr(cspyce, 'gdpool', gdpool)
    monkeypatch.setattr(cspyce, 'furnsh', furnsh)
    monkeypatch.setattr(galileo_module, 'get_spice_filecache_prefix',
                        lambda: FakeFileCachePrefix(fake.retrieved))
    monkeypatch.setattr(oops.spice, 'load_leap_seconds', load_leap_seconds)
    monkeypatch.setattr(oops.Body, 'define_solar_system', staticmethod(define_solar_system))
    monkeypatch.setattr(oops.path, 'SpicePath', spice_path)
    monkeypatch.setattr(oops.frame, 'SpiceType1Frame', spice_type1_frame)

    monkeypatch.setattr(oops.Path, '_PATH_REGISTRY', dict(oops.Path._PATH_REGISTRY))
    monkeypatch.setattr(oops.Path, '_PATH_CACHE', dict(oops.Path._PATH_CACHE))
    monkeypatch.setattr(oops.Frame, '_FRAME_REGISTRY', dict(oops.Frame._FRAME_REGISTRY))
    monkeypatch.setattr(oops.Frame, '_FRAME_CACHE', dict(oops.Frame._FRAME_CACHE))
    return fake


@pytest.fixture(autouse=True)
def reset_host() -> Any:
    """Start and end every test with the Galileo and SSI class states reset."""

    GalileoSSI._reset()
    yield
    GalileoSSI._reset()


@pytest.fixture
def ganymede_path(monkeypatch: pytest.MonkeyPatch) -> oops.Path:
    """A fixed path registered as GANYMEDE, with a Body of that name that a TrackerFrame
    can track; the Body registry is restored afterward."""

    monkeypatch.setattr(oops.Body, 'BODY_REGISTRY', dict(oops.Body.BODY_REGISTRY))
    path = oops.path.FixedPath((1.e9, 0., 0.), oops.Path.SSB, oops.Frame.J2000,
                               path_id='GANYMEDE')
    oops.Body('GANYMEDE', path, oops.Frame.J2000)
    return path


@pytest.fixture
def initialized_ssi(fake_spicedb: FakeSpiceDB) -> FakeSpiceDB:
    """GalileoSSI initialized with default options against the fake spicedb."""

    GalileoSSI.initialize()
    return fake_spicedb


def write_text(path: pathlib.Path, lines: list[str]) -> pathlib.Path:
    """Write PDS3 label lines with CRLF line ends; return the path."""

    path.write_text('\r\n'.join([*lines, 'END', '']))
    return path

##########################################################################################
