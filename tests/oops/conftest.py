##########################################################################################
# tests/oops/conftest.py
##########################################################################################
"""Shared fixtures that let the Galileo host run without SPICE kernels or spicedb.

Every spicedb function the package calls is replaced by :class:`FakeSpiceDB`, and the
oops and cspyce calls that would read SPICE kernels are replaced by kernel-free
equivalents. Everything else, including oops itself, runs for real.

The writers at the end build small synthetic SSI products, indexes, and supplements.
"""

import pathlib
import re
from dataclasses import dataclass, field
from typing import Any, cast

import cspyce
import julian
import numpy as np
import oops
import pytest
import spicedb
import vicar

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
    `cspyce.furnsh` only records the path. The Path and Frame registries, caches, and
    waypoint and wayframe pools are copied so that paths and frames made by one test never
    reach another.
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
        return oops.frame.Cmatrix(oops.Matrix3.IDENTITY, oops.Frame.J2000,
                                  frame_id=frame_id)

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
    monkeypatch.setattr(oops.Body, 'define_solar_system',
                        staticmethod(define_solar_system))
    monkeypatch.setattr(oops.path, 'SpicePath', spice_path)
    monkeypatch.setattr(oops.frame, 'SpiceType1Frame', spice_type1_frame)

    monkeypatch.setattr(oops.Path, '_PATH_REGISTRY', dict(oops.Path._PATH_REGISTRY))
    monkeypatch.setattr(oops.Path, '_PATH_CACHE', dict(oops.Path._PATH_CACHE))
    monkeypatch.setattr(oops.Frame, '_FRAME_REGISTRY', dict(oops.Frame._FRAME_REGISTRY))
    monkeypatch.setattr(oops.Frame, '_FRAME_CACHE', dict(oops.Frame._FRAME_CACHE))

    # A path or frame defined as an earlier one was shares its waypoint or wayframe, and
    # so its ID, so each class's pool is copied too; otherwise a frame that one test makes
    # could come back as another test's frame of another name. This is New Horizons'
    # fixture, kept until oops offers one helper (SETI/rms-oops#250). The pools are
    # WeakValueDictionary objects, from which garbage collection can remove a member
    # while dict(pool) looks up its keys (KeyError); items() iterates under the pool's own
    # guard, which holds such removals back.
    for path_class in oops.Path._PATH_SUBCLASSES:
        if '_WAYPOINTS' in vars(path_class):
            monkeypatch.setattr(path_class, '_WAYPOINTS',
                                dict(path_class._WAYPOINTS.items()))
    for frame_class in oops.Frame._FRAME_SUBCLASSES:
        for pool in ('_WAYFRAMES', '_FOR_NAME', '_FRAME_LOOKUP'):
            if pool in vars(frame_class):
                monkeypatch.setattr(frame_class, pool,
                                    dict(getattr(frame_class, pool).items()))
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
# An SSI product: its PDS3 label and VICAR image file
##########################################################################################

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


# The VICAR header keywords of C0349632100R that the host reads
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


def write_vicar_image(directory: pathlib.Path, name: str = 'C0349632100R.IMG', /, *,
                      array: np.ndarray | None = None,
                      **changes: float | str | list[int] | None) -> pathlib.Path:
    """Write a VICAR image file; return its path.

    The array is `array`, or :func:`image_array` if it is None. The header holds the
    VICAR keywords of C0349632100R that the host reads, changed or extended by `changes`;
    a value of None removes the keyword.
    """

    vic = vicar.VicarImage.from_array(image_array() if array is None else array)
    for (key, value) in {**_VICAR_LABEL, **changes}.items():
        if value is not None:
            vic[key] = value
    path = directory / name
    vic.write_file(path)
    return path


def write_pds3_label(directory: pathlib.Path, *, image_name: str = 'C0349632100R.IMG',
                     record: int = 1, lines: int = SSI_PIXELS,
                     extra_lines: tuple[str, ...] = (), **changes: Any) -> pathlib.Path:
    """Write the detached PDS3 label C0349632100R.LBL; return its path.

    The label's IMAGE pointer names `image_name` at `record`, counted from one, and its
    IMAGE object has `lines` lines of 800 samples. Any `extra_lines` are added as written,
    and the label's other keywords are changed by `changes` as :func:`label_dict` changes
    them.
    """

    text = ['PDS_VERSION_ID = PDS3', 'RECORD_TYPE = FIXED_LENGTH',
            f'RECORD_BYTES = {SSI_PIXELS}', f'FILE_RECORDS = {SSI_PIXELS}',
            f'^IMAGE = ("{image_name}", {record})', *extra_lines]
    text += pds3_lines(label_dict(**changes))
    text += ['OBJECT = IMAGE', f'  LINES = {lines}', f'  LINE_SAMPLES = {SSI_PIXELS}',
             '  SAMPLE_BITS = 8', '  SAMPLE_TYPE = UNSIGNED_INTEGER',
             'END_OBJECT = IMAGE']
    return write_text(directory / 'C0349632100R.LBL', text)


def write_product(directory: pathlib.Path, *, header: dict[str, Any] | None = None,
                  **changes: Any) -> pathlib.Path:
    """Write the VICAR image C0349632100R.IMG and its PDS3 label; return the label path.

    The VICAR header is changed by `header` as :func:`write_vicar_image` changes it, and
    the label by `changes` as :func:`write_pds3_label` changes it. The label's IMAGE
    pointer names the first record after the VICAR header, as an archive label's does.
    """

    image = write_vicar_image(directory, **(header or {}))
    record = cast(int, vicar.VicarLabel(image)['LBLSIZE']) // SSI_PIXELS + 1
    return write_pds3_label(directory, **{'record': record, **changes})

##########################################################################################
# An SSI index and its supplement
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

# The index declares a zero exposure and an IMAGE_TIME of UNK as null, as the archive does
_INDEX_NULLS = {'EXPOSURE_DURATION': '0', 'IMAGE_TIME': '"UNK"'}

INDEX_ROWS = [
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

SUPPLEMENT_ROWS = [
    ('GO_0017', 'J0/OPNAV/C0346405900R.LBL', 'IM8', [-1, -1, -1, -1], 'NORMAL',
     '03464059:00:0:0'),
    ('GO_0017', 'G1/GANYMEDE/C0349632100R.LBL', 'HIS', [2, 1, 399, 392], 'NORMAL',
     '03496321:00:0:0'),
    ('GO_0017', 'G1/GANYMEDE/C0349632200R.LBL', 'HIS', [2, 1, 399, 392], 'NORMAL',
     '03496322:00:0:0'),
]

INDEX_KEYWORDS: dict[str, Any] = {'MISSION_NAME': 'GALILEO',
                                  'INSTRUMENT_HOST_NAME': 'GALILEO',
                                  'INSTRUMENT_ID': 'SSI'}


def write_table(directory: pathlib.Path, name: str, object_name: str,
                columns: list[tuple[str, str, int]], rows: list[tuple[Any, ...]],
                nulls: dict[str, str] | None = None, **keywords: Any) -> pathlib.Path:
    """Write an ASCII PDS3 table and its label; return the label's path.

    An ASCII_INTEGER column holds a list of four values, each in a field of `width`
    bytes, and declares a NULL_CONSTANT of -1; an ASCII_REAL column holds one value. A
    column named in `nulls` declares that text as its NULL_CONSTANT.
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
        if nulls and column in nulls:
            text += [f'    NULL_CONSTANT = {nulls[column]}']
        if kind == 'CHARACTER':
            text += [f'    START_BYTE = {start + 1}', f'    BYTES = {width}']
            start += width + 3
        elif kind == 'ASCII_INTEGER':
            text += [f'    START_BYTE = {start}', f'    BYTES = {4 * width + 3}',
                     '    ITEMS = 4', f'    ITEM_BYTES = {width}',
                     f'    ITEM_OFFSET = {width + 1}', '    NULL_CONSTANT = -1']
            start += 4 * (width + 1)
        else:
            text += [f'    START_BYTE = {start}', f'    BYTES = {width}']
            start += width + 1
        text += ['  END_OBJECT = COLUMN']
    text += [f'END_OBJECT = {object_name}']
    return write_text(directory / f'{name}.lbl', text)


def write_index(directory: pathlib.Path, rows: list[tuple[Any, ...]] = INDEX_ROWS,
                name: str = 'GO_0017_index', **keywords: Any) -> pathlib.Path:
    """Write an SSI index and its label; return the label's path.

    The label's top-level keywords are :data:`INDEX_KEYWORDS`, changed by `keywords`.
    """

    return write_table(directory, name, 'IMAGE_INDEX_TABLE', _INDEX_COLUMNS, rows,
                       nulls=_INDEX_NULLS, **{**INDEX_KEYWORDS, **keywords})


def write_supplement(directory: pathlib.Path,
                     rows: list[tuple[Any, ...]] = SUPPLEMENT_ROWS,
                     name: str = 'GO_0017_supplemental_index') -> pathlib.Path:
    """Write the supplemental index for :func:`write_index`; return its label's path."""

    return write_table(directory, name, 'SUPPLEMENTAL_INDEX_TABLE', _SUPPLEMENT_COLUMNS,
                       rows, **INDEX_KEYWORDS)


def write_volume(directory: pathlib.Path) -> pathlib.Path:
    """Write the SSI index and its supplemental index; return the index label's path."""

    write_supplement(directory)
    return write_index(directory)

##########################################################################################
