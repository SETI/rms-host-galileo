====================
Reading a PDS3 Index
====================

:meth:`~host_galileo.oops.ssi.GalileoSSI.from_index` makes one
:class:`~oops.observation.Snapshot` for each row of a PDS3 index. It reads only the
index tables and opens none of the images, which makes it much faster than opening a
volume's images one at a time.

It returns an iterator, which gives the Snapshots row by row, in the order of the index.
The call itself checks the arguments and reads the index, so a wrong argument or an index
that cannot be read raises an error at once. Each row's Snapshot is built only when the
iteration reaches that row, so an error in a row is raised there, after the rows before
it. The iteration ends at that row: the rows after it are not given. The ``filter``
argument, described below, can leave such a row out. To build every row at once, pass
the iterator to ``list()``:

.. code-block:: python

   from host_galileo import GalileoSSI

   index = 'metadata/GO_0xxx/GO_0017/GO_0017_index.lbl'
   for snapshot in GalileoSSI.from_index(index):
       print(snapshot.basename, snapshot.time)

   snapshots = list(GalileoSSI.from_index(index))

The index and its supplement
============================

The metadata tree of the Ring-Moon Systems Node holdings gives each GO volume a PDS3
index, such as ``GO_0017_index.lbl``, with one row per SSI image, and a supplemental
index beside it, ``GO_0017_supplemental_index.lbl``, whose rows give further values
for the same images. The host needs both. Each row of the supplement is
merged into the matching row of the index; where both tables give a column, the index's
value is kept.

A volume also carries its own PDS3 index, in its ``INDEX`` directory, such as
``GO_0017/INDEX/IMGINDEX.LBL``. That index has no supplement, so it lacks the telemetry
format and the cutout window of each image, and the host does not read it. An index in
a directory named ``INDEX``, in any case, raises :exc:`~oops.host.HostIndexError`, with a
message that points to the metadata index instead.

The ``supplement`` argument chooses the supplemental index:

- By default, the host reads the file beside the index whose name has
  ``_supplemental_index`` in place of ``_index``, in upper or lower case. That file must
  exist.
- A path or URL reads that file instead.
- None raises :exc:`~oops.host.HostIndexError`, unless ``row_dicts`` is given.

The two tables must list the same images, but not necessarily in the same order. When
they list them in the same order, as every GO volume does, each index row is merged
with the supplement row in the same position. When they do not, each index row is merged
with the supplement row that has the same ``VOLUME_ID`` and ``FILE_SPECIFICATION_NAME``.
The host raises :exc:`~oops.host.HostIndexError` if an image appears in only one of the
tables, or if, when the orders differ, either table lists an image twice. The error
lists the first few such images and both tables.

Other arguments
===============

``label``
    The index's PDS3 label, if it has already been parsed.

``row_dicts``
    The rows themselves, one dictionary per row, if the tables have already been read.
    The host then reads neither table, and ``supplement`` and ``label`` are ignored.

``filter``
    A dictionary of column names and the values to keep. A row is kept when each column
    equals its value, is one of a set of values, or falls in a range given as a tuple of
    two values. A column that no row has raises :exc:`~oops.OopsValueError`, and a
    range that cannot be compared with a column's values raises
    :exc:`~oops.OopsTypeError`, both at once.

``select``
    Must be None. Each row is one image, so there is nothing to select, and any other
    value raises :exc:`~oops.OopsValueError`.

``full_fov``, ``return_all_planets``, ``pds3_options``, and the standard overrides
    The same as for :meth:`~host_galileo.oops.ssi.GalileoSSI.from_file`, applied to every
    row.

Any other keyword raises :exc:`~oops.OopsTypeError`, so that a misspelled argument is
not silently ignored.

What one row gives
==================

Each row gives a :class:`~oops.observation.Snapshot` with no ``data``. Its ``dict``
subfield holds the merged row, with the same treatment of missing values as a label (see
:doc:`user_guide_metadata`). The placeholders that the tables write for a value an image
does not have are left out: the supplement's ``CUT_OUT_WINDOW`` and ``TRUTH_WINDOW`` of
-1, -1, -1, -1, its ``ICT_DESPIKE_THRESHOLD`` and ``ICT_QUANTIZATION_STEP_SIZE`` of -1,
a compression ratio of 0, the ``SMEAR_MAGNITUDE`` of 0.1 that the index gives every
image of the Jupiter tour, and the ``ORBIT_NUMBER`` of -99 that the index gives an image
taken before the tour. An exposure of 0 is a real value and is kept.

The indexes of GO_0017 through GO_0023, the only ones that give compression ratios,
write each image's three ratios in the wrong columns. The value in
``ENCODING_COMPRESSION_RATIO`` is the image's maximum ratio, the value in
``ENCODING_MIN_COMPRESSION_RATIO`` is its ratio, and the value in
``ENCODING_MAX_COMPRESSION_RATIO`` is its minimum. For example, the product label of
``GO_0021/C21/IO/C0506406074R`` gives a ratio of 3.758, a minimum of 3.03, and a maximum
of 5.428, while the index gives 5.428, 3.758, and 3.03. The host gives each value under
the name of the ratio it is, as the product label does.

The same indexes give every image of the Jupiter tour the data set version 1.0, as in
``GO-J/JSA-SSI-2-REDR-V1.0``, while every product label gives version 1.1. The host
gives ``DATA_SET_ID`` as the product label does.

The subfields are those of an image read from its label (see
:doc:`user_guide_metadata`), except ``data``, ``label_dict``, and ``header_dict``. The
telemetry format comes from the supplement, so a row passed in ``row_dicts`` without the
supplement's columns gets the full field of view.

Every row gives a Snapshot, including a row whose exposure is zero. Building a row's
Snapshot raises no error about pointing. If the C kernels hold no pointing for an image,
its Snapshot's frame raises :exc:`oops.OopsSpiceError` whenever it is evaluated, with
the same message that :meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` gives (see
:doc:`user_guide_errors`).

File paths
==========

A row's ``FILE_SPECIFICATION_NAME`` is the path of the PDS3 label, but each Snapshot
refers to the image file, with ``.IMG`` in place of ``.LBL`` in the same case.
``filepath``, ``basename``, and ``image_url`` give it as a full path or URL, formed from
the location of the index:

- For an index in the ``metadata`` tree of the Ring-Moon Systems Node holdings, such as
  ``.../metadata/GO_0xxx/GO_0017/GO_0017_index.lbl``, the volume directories are in the
  matching ``volumes`` tree, ``.../volumes/GO_0xxx``.
- For an index anywhere else, they are in the index's own directory.

The row's ``VOLUME_ID`` and file specification complete the path. ``abspath`` is present
only when the image file is local or already in the local cache. Reading an index never
creates a directory, even when the volumes are not there.

Reading an index through oops.Host
==================================

:meth:`oops.host.Host.from_index` recognizes an SSI index from its label and passes it,
with every other argument, to :meth:`~host_galileo.oops.ssi.GalileoSSI.from_index`,
which reads and pairs the two tables itself. The two give the same Snapshots. A
volume's own index does not identify its instrument in its label, so
:meth:`oops.host.Host.from_index` recognizes it from its first row and passes it to the
host, which raises :exc:`~oops.host.HostIndexError` as described above.
