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
from collections.abc import Iterator
from typing import Any

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

class GalileoSSI(Host):
    # The runtime names the first parameter `fileinfo`, where Host names it `filepath`.
    # Given the image file, only it is read; given a PDS3 label, the label and its image
    # file are read and their metadata merged, except that with astrometry=True the label
    # is read alone and the image file need not exist. A value that no label read gives
    # has no subfield and no dict key; each item that both give takes its value from one
    # fixed label. A label or VICAR header that cannot be parsed, or a label with no
    # IMAGE pointer, raises oops.host.PdsHostError or VicarHostError (subclasses of
    # HostFormatError), as does an image for which no label read gives the exposure, the
    # time, or the telemetry format; a missing image file behind a label raises
    # FileNotFoundError. Unless a frame is given, an image whose pointing the C kernels do
    # not hold raises oops.OopsSpiceError, whose cause is the SPICE exception.
    @staticmethod
    def from_file(  # type: ignore[override]
        fileinfo: str | pathlib.Path | FCPath | Pds3Label | VicarLabel
        | VicarImage, *,
        full_fov: bool = False, return_all_planets: bool = False,
        pds3_options: dict[str, Any] | None = None, astrometry: bool = False,
        target: str | Body | None = None, lightsource: LightSource | None = None,
        path: str | Path | None = None, frame: str | Frame | None = None,
        cadence: Cadence | None = None, fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None, **kwargs: Any) -> Snapshot: ...
    # Every row is one Snapshot. The supplement is required: None raises
    # oops.host.HostIndexError unless row_dicts is given. select must be None
    # (oops.OopsValueError), and any other keyword raises oops.OopsTypeError. The
    # iterator builds each row's Snapshot when it reaches the row.
    @staticmethod
    def from_index(
        filepath: str | pathlib.Path | FCPath, *,
        supplement: str | pathlib.Path | FCPath | object | None = ...,
        label: Pds3Label | None = None,
        row_dicts: list[dict[str, Any]] | None = None,
        pds3_options: dict[str, Any] | None = None,
        filter: dict[str, Any] | None = None, full_fov: bool = False,
        return_all_planets: bool = False,
        select: int | slice | str | tuple[int | slice | str, ...] | None = None,
        target: str | Body | None = None, lightsource: LightSource | None = None,
        path: str | Path | None = None, frame: str | Frame | None = None,
        cadence: Cadence | None = None, fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None, **kwargs: Any) -> Iterator[Snapshot]: ...
    @staticmethod
    def initialize(*, planets: list[Any] | None = None, asof: str | None = None,
        irregulars: bool = True) -> None: ...

##########################################################################################
