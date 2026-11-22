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
The Ring-Moon Systems Node publishes these resources, the kernels and ``SPICE.db``
included, at ``https://storage.googleapis.com/rms-node-oops-resources``, a URL that
``OOPS_RESOURCES`` can name directly.

The host furnishes one fixed set of kernels, which covers the whole mission. They lie
under ``General`` and ``Galileo`` in the ``SPICE_PATH`` directory: the leap-second
kernel, the SSI instrument kernel, the frame kernel, the spacecraft clock kernel, the
planetary and spacecraft ephemerides, and the scan platform C kernels. Other kernels
are chosen through the database before the fixed set is furnished: a leap-second
kernel; the SSI instrument kernel, with the frame and spacecraft clock kernels that
belong to it; and the kernels that OOPS loads to define the planets and moons, which
are the planetary constants and the ephemerides of the planets and their moons. Where
two kernels give the same information, SPICE uses the one furnished last, so wherever
the fixed set and the database's choices overlap, the fixed set is what counts.

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
camera's fields of view. To set the options described below, call
:meth:`~host_galileo.oops.ssi.GalileoSSI.initialize` before opening the first image.
Later calls to :meth:`~host_galileo.oops.ssi.GalileoSSI.initialize` are ignored.

The ``planets`` option chooses the bodies that OOPS defines. It takes one planet
number, from 1 for Mercury to 9 for Pluto, or a sequence of them. This call defines
only Jupiter, planet 5, and its moons:

.. code-block:: python

   from host_galileo import GalileoSSI

   GalileoSSI.initialize(planets=(5,))

The ``asof`` option takes a date, such as ``'1997-01-01'``, and each kernel is then
chosen from those released no later than the start (00:00) of that date. Where every
kernel that a query could choose was released after the start of the date, the query
silently uses the kernel it would choose without a date. That happens only when the
whole query comes back empty: a query that finds some kernels released by the start of
the date uses just those. The date changes less than that suggests:

- The fixed set of kernels is loaded without the date, and last, so it overrides the
  database's choices wherever the two overlap. The SSI instrument, frame, and
  spacecraft clock kernels, the C kernels, and every ephemeris in the fixed set are the
  same whatever the date. For example, the positions of Jupiter, Io, and Europa seen
  from Galileo are the same with ``asof='1997-01-01'`` as without a date.
- The leap seconds come from the fixed set's leap seconds kernel, naif0012, which is
  loaded last, so they are the same whatever the date.

The date can change the planetary constants and the ephemerides of bodies that the
fixed set does not cover. It could also change the camera's fields of view, which are
built from the instrument kernel that the database chooses, before the fixed set is
furnished. The database lists a single SSI instrument kernel, the one in the fixed
set, so the fields of view are the same whatever the date.

Image files and indexes
=======================

Each image file, PDS3 label, or index is passed by its path or URL, so the location of
the archive is given with each call rather than in an environment variable. A remote
file is retrieved and cached when it is read.
