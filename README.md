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
observation geometry. It reads images from Galileo's Solid State Imaging (SSI) camera, as
the Planetary Data System (PDS) archives them, and turns each one into an OOPS `Snapshot`:
the image array, together with the time, field of view, pointing, target, and light
source that OOPS needs to compute the image's geometry. It is a product of the
[PDS Ring-Moon Systems Node](https://pds-rings.seti.org).

# Features

- **Either file of a product**: an image can be opened from its VICAR image file or from
  its PDS3 label, and the pixels are identical either way.
- **Metadata from both labels**: the values in the PDS3 label and in the VICAR header are
  merged under the PDS3 names, with one fixed choice for each value on which they
  disagree.
- **Whole volumes at once**: a PDS3 index and its supplemental index give one `Snapshot`
  per image, without opening the images.
- **Pointing from SPICE**: the pointing comes from the Galileo C kernels, or from a frame
  that the caller supplies.
- **Exposure timing**: an image that has no image time is timed from its spacecraft clock
  count and the camera's shutter delay.

# Installation

`rms-host-galileo` needs Python 3.11 or later. Install it from PyPI with:

```sh
pip install rms-host-galileo
```

The host also needs the Galileo SPICE kernels. The
[setup chapter](https://rms-host-galileo.readthedocs.io/en/latest/user_guide/user_guide_setup.html)
of the user guide describes the kernels and the environment variables that locate them.

# Quick Start

Importing `host_galileo.oops` registers the host with `oops.Host`, which can then open
any SSI image:

```python
import oops
import host_galileo.oops

obs = oops.Host.from_file('GO_0017/G1/GANYMEDE/C0349632100R.IMG')
```

The host class can also be used directly, for one image or for every image in a volume:

```python
from host_galileo import GalileoSSI

snapshot = GalileoSSI.from_file('GO_0017/G1/GANYMEDE/C0349632100R.LBL')
snapshots = GalileoSSI.from_index('GO_0017_index.lbl')
```

Replace each path with the location of the file in your copy of the archive, or with
its URL.

# Documentation

The [documentation](https://rms-host-galileo.readthedocs.io) holds the user guide and
the API reference. To build it locally, run `make html` in the `docs` directory after
installing the package with its `docs` extra.

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
