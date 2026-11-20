=====================
Reading an Image File
=====================

:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` makes one
:class:`~oops.observation.Snapshot` from one SSI image. It accepts either file of the
product, the VICAR image file or its PDS3 label:

.. code-block:: python

   from host_galileo import GalileoSSI

   snapshot = GalileoSSI.from_file('GO_0017/G1/GANYMEDE/C0349632100R.IMG')
   snapshot = GalileoSSI.from_file('GO_0017/G1/GANYMEDE/C0349632100R.LBL')

What is read
============

What the host reads depends on the file it is given and on the ``astrometry`` option.
An astrometry read is for work that needs only the image's time and geometry, so it
skips the image array.

.. list-table::
   :header-rows: 1

   * - Input
     - Normal read
     - With ``astrometry=True``
   * - VICAR image file
     - The image file alone: its header and its array.
     - The image file's header alone.
   * - PDS3 label
     - The label, and the whole image file beside it.
     - The label alone. The image file need not exist.

The metadata come from whatever was read, as :doc:`user_guide_metadata` describes. The
``data`` subfield holds the image array except in an astrometry read.

The host always finds the image array through the VICAR file's own header, never
through the record number in the label's ``^IMAGE`` pointer. Some labels give that
number wrongly, so this keeps the array identical from either input.

Finding the image file
======================

Given a PDS3 label, the host looks in the label's directory for the image file that the
label's ``^IMAGE`` pointer gives. It tries the name as the label spells it, then
in the other case, so that a copy of the archive with lower-case file names still works.
An astrometry read does not look for the file, so its path subfields give the name as
the label spells it.
A path that ends in ``.lbl``, in any case, is read as a PDS3 label.

The ``filepath``, ``basename``, ``abspath``, and ``image_url`` subfields always give the
image file, also when the input is the label. ``abspath`` is the image file's path on
the local disk, so it is present only when the file is local or is already in the local
cache of a remote file. After an astrometry read of a remote label, the image file has
not been retrieved, so ``abspath`` is absent unless an earlier read put the file in the
cache. Opening an image leaves the directories of the file system as they are.

Parsed input
============

A label or image that has already been parsed can be passed in its place: a
:class:`~pdsparser.Pds3Label` is read as a PDS3 label, and a :class:`~vicar.VicarLabel`
or :class:`~vicar.VicarImage` as the image file.

Parsed content that was read from text, rather than from a file, has no file path. Its
Snapshot has no ``filepath``, ``basename``, ``abspath``, or ``image_url`` subfield. An
astrometry read of it works as usual, and so does a normal read of a
:class:`~vicar.VicarImage`, which already holds its array. A normal read of a
:class:`~pdsparser.Pds3Label` or :class:`~vicar.VicarLabel` needs the image file, so it
raises :exc:`~oops.host.HostFormatError`, whose message says that no file path is
known.

The ``header_dict`` subfield holds the VICAR header as the file writes it. A
:class:`~vicar.VicarImage` sets its ``HOST`` and ``BHOST`` to the machine that read it,
so when a :class:`~vicar.VicarImage` that has a file path is passed in, the host reads
the header again from that file. A :class:`~vicar.VicarImage` without a file path, such
as one made by ``VicarImage.from_file()``, which keeps no path, has only the header that
it holds, so its ``header_dict`` gives the reading machine as ``HOST`` and ``BHOST``.

Parsed input can also be passed to :meth:`oops.host.Host.from_file`, which identifies it
as an SSI label or image and passes it on to
:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` as it is.

Options
=======

These options are particular to the Galileo host:

``full_fov``
    False, the default, trims the image and the field of view to the cutout window, as
    :doc:`user_guide_products` describes. True keeps the full frame.

``return_all_planets``
    True lists in the ``spice_kernels`` subfield the kernels of every planet, and not
    only those of the planet that the mission was observing at the time.

``pds3_options``
    Options for the PDS3 label parser. The host starts with the parser's fast method and,
    if that cannot read a label, tries the strict method and then the loose one. A method
    given in ``pds3_options`` starts that sequence at the method given.

``astrometry``
    True for an astrometry read, as described above.

``select``
    Must be None. An SSI image is one observation, so there is nothing to select, and
    any other value raises :exc:`~oops.OopsValueError`.

The standard OOPS overrides, ``target``, ``lightsource``, ``path``, ``frame``,
``cadence``, ``fov``, ``calibrations``, ``timeshift``, ``navigation``, ``tracker``, and
``parallel``, replace what the host would otherwise build. They are listed with
:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file`. Any other keyword raises
:exc:`~oops.OopsTypeError`, so that a misspelled argument is caught.

Reading a label through oops.Host
=================================

:meth:`oops.host.Host.from_file` only identifies the host. It parses a PDS3 label with
the parser's fast method, then the strict one, then the loose one, until one of them
reads it, and it reads a VICAR file's header alone. It then passes the original path,
with the caller's options and overrides, to
:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file`. Any SSI label or image file that
the host reads directly therefore also reads through :meth:`oops.host.Host.from_file`,
and the two give the same Snapshot.
