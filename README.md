<!-- pyml disable MD025 -->
# rms-host-galileo

[![GitHub release; latest by date](https://img.shields.io/github/v/release/SETI/rms-host-galileo)](https://github.com/SETI/rms-host-galileo/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/SETI/rms-host-galileo)](https://github.com/SETI/rms-host-galileo/releases)
[![Test Status](https://img.shields.io/github/actions/workflow/status/SETI/rms-host-galileo/run-tests.yml?branch=main)](https://github.com/SETI/rms-host-galileo/actions)
[![Documentation Status](https://readthedocs.org/projects/rms-host-galileo/badge/?version=latest)](https://rms-host-galileo.readthedocs.io/en/latest/?badge=latest)
[![Code coverage](https://img.shields.io/codecov/c/github/SETI/rms-host-galileo/main?logo=codecov)](https://codecov.io/gh/SETI/rms-host-galileo)
<br />
[![PyPI - Version](https://img.shields.io/pypi/v/rms-host-galileo)](https://pypi.org/project/rms-host-galileo)
[![PyPI - Format](https://img.shields.io/pypi/format/rms-host-galileo)](https://pypi.org/project/rms-host-galileo)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/rms-host-galileo)](https://pypi.org/project/rms-host-galileo)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/rms-host-galileo)](https://pypi.org/project/rms-host-galileo)
<br />
[![GitHub commits since latest release](https://img.shields.io/github/commits-since/SETI/rms-host-galileo/latest)](https://github.com/SETI/rms-host-galileo/commits/main/)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/m/SETI/rms-host-galileo)](https://github.com/SETI/rms-host-galileo/commits/main/)
[![GitHub last commit](https://img.shields.io/github/last-commit/SETI/rms-host-galileo)](https://github.com/SETI/rms-host-galileo/commits/main/)
<br />
[![Number of GitHub open issues](https://img.shields.io/github/issues-raw/SETI/rms-host-galileo)](https://github.com/SETI/rms-host-galileo/issues)
[![Number of GitHub closed issues](https://img.shields.io/github/issues-closed-raw/SETI/rms-host-galileo)](https://github.com/SETI/rms-host-galileo/issues)
[![Number of GitHub open pull requests](https://img.shields.io/github/issues-pr-raw/SETI/rms-host-galileo)](https://github.com/SETI/rms-host-galileo/pulls)
[![Number of GitHub closed pull requests](https://img.shields.io/github/issues-pr-closed-raw/SETI/rms-host-galileo)](https://github.com/SETI/rms-host-galileo/pulls)
<br />
![GitHub License](https://img.shields.io/github/license/SETI/rms-host-galileo)
[![Number of GitHub stars](https://img.shields.io/github/stars/SETI/rms-host-galileo)](https://github.com/SETI/rms-host-galileo/stargazers)
![GitHub forks](https://img.shields.io/github/forks/SETI/rms-host-galileo)
<!-- start-after-point -->

# Introduction

`rms-host-galileo` is the Galileo mission host for
[OOPS](https://github.com/SETI/rms-oops), the Object-Oriented Python SPICE toolkit for
observation geometry. It turns Galileo Solid State Imaging (SSI) images into OOPS
`Snapshot` objects, with the frame, field of view, timing, target, and light source that
each image requires. It is a product of the
[PDS Ring-Moon Systems Node](https://pds-rings.seti.org).

# Installation

The `host_galileo` module is available via the `rms-host-galileo` package on PyPI and can
be installed with:

```sh
pip install rms-host-galileo
```

The host needs the `oops.Host` interface, which is not yet in a released version of
`rms-oops`. Until it is, install `rms-oops` from its development branch:

```sh
pip install "rms-oops @ git+https://github.com/SETI/rms-oops@mrs_260925_host_reorg"
```

Geometry calculations need the Galileo SPICE kernels. The host reads a fixed set of them
through the SPICE file cache of `spicedb`, located by the `SPICE_PATH` environment
variable (or `$OOPS_RESOURCES/SPICE`), and the SSI instrument kernel through the kernel
database named by `SPICE_SQLITE_DB_NAME`.

# Getting Started

Importing `host_galileo.oops` registers the Galileo host with `oops.Host`, which can then
identify and read any SSI image, from either its PDS3 label or its VICAR file:

```python
import oops
import host_galileo.oops

obs = oops.Host.from_file('GO_0017/G1/GANYMEDE/C0349632100R.LBL')
```

The host class can also be used directly. It reads a volume index together with the
supplemental index beside it (here `GO_0017_supplemental_index.lbl`), which must exist:

```python
from host_galileo import GalileoSSI

snapshot = GalileoSSI.from_file('GO_0017/G1/GANYMEDE/C0349632100R.IMG')
snapshots = GalileoSSI.from_index('GO_0017_index.lbl')
```

# Features

- `GalileoSSI.from_file()` makes a `Snapshot` from an SSI image. The image and field of
  view are trimmed to the cutout window given in the label unless `full_fov=True`, and
  `astrometry=True` reads no image data.
- `GalileoSSI.from_index()` makes one `Snapshot` for each row of a volume index,
  including rows with a zero exposure. It merges in the columns of the supplemental
  index, matching its rows to the index's by VOLUME_ID and FILE_SPECIFICATION_NAME.
- The exposure is centered on the label's IMAGE_TIME. When the label gives none, the
  time comes from the spacecraft clock count and the shutter delay for the frame
  duration.
- A target named in the label is repaired where the label is known to be wrong, and an
  image of a star uses that star as its light source.

Details are available in the
[module documentation](https://rms-host-galileo.readthedocs.io/en/latest/module.html).

# Contributing

Information on contributing to this package can be found in the
[Contributing Guide](https://github.com/SETI/rms-host-galileo/blob/main/CONTRIBUTING.md).

# Links

- [Documentation](https://rms-host-galileo.readthedocs.io)
- [Repository](https://github.com/SETI/rms-host-galileo)
- [Issue tracker](https://github.com/SETI/rms-host-galileo/issues)
- [PyPi](https://pypi.org/project/rms-host-galileo)

# Licensing

This code is licensed under the [Apache License v2.0](https://github.com/SETI/rms-host-galileo/blob/main/LICENSE).
