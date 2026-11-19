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
reads no image array.

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
A path that ends in ``.lbl``, in any case, is read as a PDS3 label.

The ``filepath``, ``basename``, ``abspath``, and ``image_url`` subfields always give the
image file, also when the input is the label. After an astrometry read of a remote
label, the image file has not been retrieved, so ``abspath`` is absent unless an earlier
read put the file in the local cache.

Parsed input
============

A label or image that has already been parsed can be passed in its place: a
:class:`~pdsparser.Pds3Label` is read as a PDS3 label, and a :class:`~vicar.VicarLabel`
or :class:`~vicar.VicarImage` as the image file.

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

The standard OOPS overrides, ``target``, ``lightsource``, ``path``, ``frame``,
``cadence``, ``fov``, ``calibrations``, ``timeshift``, ``navigation``, ``tracker``, and
``parallel``, replace what the host would otherwise build. They are listed with
:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file`. Any other keyword is accepted and
ignored.
