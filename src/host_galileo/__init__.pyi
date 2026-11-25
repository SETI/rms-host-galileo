##########################################################################################
# host_galileo/__init__.pyi
##########################################################################################
"""Type stub for :mod:`host_galileo`.

The source types its methods in their docstrings rather than in their signatures, so the
type information for public symbols is published here instead. Only package stubs exist,
so a name is annotated when it is imported from the package that exports it and not when
it is imported from the module that defines it. The stub describes the shape of the API
exactly: every public name, its parameters, which of them are keyword-only, and which have
defaults. Types are given where they are unambiguous and are `Any` elsewhere.

The detectors and support methods that a Host subclass overrides or calls, such as
`_register` and `_detect_in_pds3`, are private and are not published here.
"""

import pathlib
from collections.abc import Iterable, Iterator, Mapping
from typing import Any, overload

from filecache import FCPath
from pdsparser import Pds3Label
from vicar import VicarImage, VicarLabel
from oops import Body
from oops.cadence import Cadence
from oops.calibration import Calibration
from oops.fov import FOV
from oops.frame import Frame
from oops.host import Host
from oops.lightsource import LightSource
from oops.observation import Observation, Snapshot
from oops.path import Path

__all__ = ['GalileoSSI']

_Select = int | slice | str | tuple[int | slice | str, ...]

class GalileoSSI(Host):
    # The runtime names the first parameter `fileinfo`, where Host names it `filepath`,
    # and takes no `**kwargs`. Without `select`, the result is the image's one Snapshot.
    @overload  # type: ignore[override]
    @staticmethod
    def from_file(
        fileinfo: str | pathlib.Path | FCPath | Pds3Label | VicarLabel
        | VicarImage, *,
        full_fov: bool = False, return_all_planets: bool = False,
        pds3_options: dict[str, Any] | None = None, astrometry: bool = False,
        select: None = None, target: str | Body | None = None,
        lightsource: LightSource | None = None, path: str | Path | None = None,
        frame: str | Frame | None = None, cadence: Cadence | None = None,
        fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None) -> Snapshot: ...
    @overload
    @staticmethod
    def from_file(
        fileinfo: str | pathlib.Path | FCPath | Pds3Label | VicarLabel
        | VicarImage, *,
        full_fov: bool = False, return_all_planets: bool = False,
        pds3_options: dict[str, Any] | None = None, astrometry: bool = False,
        select: _Select, target: str | Body | None = None,
        lightsource: LightSource | None = None, path: str | Path | None = None,
        frame: str | Frame | None = None, cadence: Cadence | None = None,
        fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None) -> Snapshot | tuple[Snapshot, ...]: ...
    # The runtime takes no `**kwargs`, and every row is one Snapshot, so the iterator is
    # typed more narrowly than Host's without `select`.
    @overload  # type: ignore[override]
    @staticmethod
    def from_index(
        filepath: str | pathlib.Path | FCPath, *,
        supplement: str | pathlib.Path | FCPath | object | None = ...,
        label: Pds3Label | None = None,
        row_dicts: Iterable[dict[str, Any]] | None = None,
        pds3_options: dict[str, Any] | None = None,
        filter: dict[str, Any] | None = None, full_fov: bool = False,
        return_all_planets: bool = False, select: None = None,
        target: str | Body | None = None, lightsource: LightSource | None = None,
        path: str | Path | None = None, frame: str | Frame | None = None,
        cadence: Cadence | None = None, fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None) -> Iterator[Snapshot]: ...
    @overload
    @staticmethod
    def from_index(
        filepath: str | pathlib.Path | FCPath, *,
        supplement: str | pathlib.Path | FCPath | object | None = ...,
        label: Pds3Label | None = None,
        row_dicts: Iterable[dict[str, Any]] | None = None,
        pds3_options: dict[str, Any] | None = None,
        filter: dict[str, Any] | None = None, full_fov: bool = False,
        return_all_planets: bool = False, select: _Select,
        target: str | Body | None = None, lightsource: LightSource | None = None,
        path: str | Path | None = None, frame: str | Frame | None = None,
        cadence: Cadence | None = None, fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None
    ) -> Iterator[Snapshot | tuple[Snapshot, ...]]: ...
    @staticmethod
    def tstart_from_sclk(label_dict: Mapping[str, Any]) -> float: ...
    @staticmethod
    def initialize(*, planets: list[Any] | None = None, asof: str | None = None,
        irregulars: bool = True) -> None: ...

##########################################################################################
