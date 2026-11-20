=====
Setup
=====

Installation
============

The package needs Python 3.11 or later and installs from PyPI:

.. code-block:: sh

   pip install rms-host-galileo

It brings in OOPS (``rms-oops``) and the readers that the host uses for PDS3 labels,
VICAR files, and PDS3 index tables.

SPICE kernels
=============

The geometry of every image comes from the Galileo SPICE kernels. The host reads them
through ``spicedb``, which is distributed with OOPS and finds them with two environment
variables:

``SPICE_PATH``
    The directory that holds the kernels, as a local path or a URL. If it is not set,
    ``$OOPS_RESOURCES/SPICE`` is used.

``SPICE_SQLITE_DB_NAME``
    The SQLite database that lists the kernels. If it is not set, ``SPICE.db`` in the
    ``SPICE_PATH`` directory is used.

The kernels and the database are the ones that OOPS itself uses. The
`rms-oops README <https://github.com/SETI/rms-oops#environment-variables>`_ describes
these variables and the ``SPICE`` directory of ``$OOPS_RESOURCES`` that holds them.

The host furnishes one fixed set of kernels, which covers the whole mission. They lie
under ``General`` and ``Galileo`` in the ``SPICE_PATH`` directory: the leap-second
kernel, the spacecraft clock kernel, the frame kernel, the planetary and spacecraft
ephemerides, and the scan platform C kernels. The SSI instrument kernel is chosen
through the database.

Registering the host
====================

The host classes import lazily. ``import host_galileo`` gives access to
:class:`~host_galileo.oops.ssi.GalileoSSI` without importing OOPS until the class is
first used. Importing ``host_galileo.oops`` registers the host with
:class:`~oops.host.Host`, so that :meth:`oops.host.Host.from_file` recognizes an SSI
image:

.. code-block:: python

   import oops
   import host_galileo.oops

   obs = oops.Host.from_file('GO_0017/G1/GANYMEDE/C0349632100R.IMG')

Initialization
==============

The first call to :meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` or
:meth:`~host_galileo.oops.ssi.GalileoSSI.from_index` loads the kernels and defines the
camera's fields of view. To choose the bodies that OOPS defines, or to use only the
instrument kernel that existed before a given date, call
:meth:`~host_galileo.oops.ssi.GalileoSSI.initialize` before opening the first image.
Its ``planets`` option takes one planet number, from 1 for Mercury to 9 for Pluto, or a
sequence of them. This call defines Jupiter, planet 5, and its moons, and no other
planet:

.. code-block:: python

   from host_galileo import GalileoSSI

   GalileoSSI.initialize(planets=(5,))

Later calls to :meth:`~host_galileo.oops.ssi.GalileoSSI.initialize` have no effect.

Image files and indexes
=======================

The host needs no environment variable for the archive itself. Each image file, PDS3
label, or index is passed by its path or URL, and a remote file is retrieved and cached
when it is read.
