====================
Reading a PDS3 Index
====================

:meth:`~host_galileo.oops.ssi.GalileoSSI.from_index` makes one
:class:`~oops.observation.Snapshot` for each row of a PDS3 index. It reads the index
tables alone, without the images, which makes it much faster than opening a volume's
images one at a time.

It returns an iterator, which gives the Snapshots row by row, in the order of the index.
The call itself checks the arguments and reads the PDS3 labels of the index and its
supplement, so a wrong argument or a label that cannot be read raises an error at once.
The rows of the two tables are read as the iteration goes on, a few at a time, and each
row's Snapshot is built when the iteration reaches that row. The memory used stays about
the same however long the index is. An error in a row is raised when the iteration
reaches it, after the rows before it, and the iteration ends there: the rows after it
are not given. The one exception is a row longer than its label says, whose error comes
a few rows early, as described below. The ``filter`` argument, described below, can leave such a row out. To build
every row at once, pass the iterator to ``list()``:

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
for the same images. The host needs both. Each row of the supplement is merged into the
matching row of the index; where both tables give a column, the index's value is kept.

A volume also carries its own PDS3 index, in its ``INDEX`` directory, such as
``GO_0017/INDEX/IMGINDEX.LBL``. That index has no supplement, so it lacks the telemetry
format and the cutout window of each image, and the host does not read it. An index in a
directory named ``INDEX``, in upper or lower case, raises
:exc:`~oops.host.HostIndexError` before anything else is read, with a message that points
to the metadata index instead.

The ``supplement`` argument chooses the supplemental index:

- By default, the host reads the file beside the index whose name has
  ``_supplemental_index`` in place of ``_index``, in upper or lower case. That file must
  exist.
- A path or URL reads that file instead.
- None raises :exc:`~oops.host.HostIndexError`, unless ``row_dicts`` is given.

Every SSI index lists its images in the same order as its supplement, so the host reads
the two tables in step and merges each index row with the supplement row in the same
position. The two rows must give the same ``VOLUME_ID`` and ``FILE_SPECIFICATION_NAME``.
If they do not, or if one table ends before the other, the host raises
:exc:`~oops.host.HostIndexError` when the iteration reaches that row, and the message
gives the row and both tables. If the supplement has rows left over, the error comes
after the last row of the index.

A table file that does not match its label is read as follows:

- A table file that holds fewer rows than its label says gives every row it holds, then
  raises :exc:`~oops.host.HostIndexError`. One that holds more raises the same error
  after the last row that its label lists.
- A row longer than the label says raises :exc:`~oops.host.HostIndexError` at the first
  of the ten rows read with it, and none of those ten rows is given. A row shorter than
  the label says is read as it is.
- A value that cannot be read as the number its column describes is not an error of the
  table: the table reader issues a warning, and the row is given with that value as
  text. Building the row's Snapshot may then fail, for example with a
  :exc:`TypeError` when the value is the exposure duration.

Other arguments
===============

``label``
    The index's PDS3 label, if it has already been parsed.

``row_dicts``
    The rows themselves, one dictionary per row, if the tables have already been read,
    each row already merged with its supplement row. Any iterable works, and the host
    takes the rows one at a time as the iteration needs them. The host then reads
    neither table, and ``supplement`` is not used. It still reads the index's PDS3
    label, unless ``label`` is given, because the label tells whether the table writes
    the compression ratios in rotated columns (see below).

``filter``
    A dictionary of column names and the values to keep. A row is kept when each column
    equals its value, is one of a set of values, or falls in a range given as a tuple of
    two values. The columns of both tables can be used. A column that neither table has
    raises :exc:`~oops.OopsValueError` at once. With ``row_dicts``, the columns are
    checked against the first row: at once if ``row_dicts`` is a sequence such as a
    list, and when the iteration reads the first row if it is an iterator. A range that
    cannot be compared with a row's value raises :exc:`~oops.OopsTypeError` when the
    iteration reaches that row.

``select``
    Must be None. Each row is one image, so there is nothing to select, and any other
    value raises :exc:`~oops.OopsValueError`.

``full_fov``, ``return_all_planets``, ``pds3_options``, and the standard overrides
    The same as for :meth:`~host_galileo.oops.ssi.GalileoSSI.from_file`, applied to every
    row.

Any other keyword raises :exc:`~oops.OopsTypeError`, so that a misspelled argument is
caught.

What one row gives
==================

Each row gives a :class:`~oops.observation.Snapshot` with no ``data``. Its ``dict``
subfield holds the merged row, with the same treatment of missing values as a label (see
:doc:`user_guide_metadata`). A table's PDS3 label can declare, for each column, a value
that stands for "no value", and such values are left out: in the current tables, the
supplement's ``CUT_OUT_WINDOW`` and ``TRUTH_WINDOW`` of -1, -1, -1, -1, its
``ICT_DESPIKE_THRESHOLD`` and ``ICT_QUANTIZATION_STEP_SIZE`` of -1, a compression ratio
of 0, a number of -999, and the ``ORBIT_NUMBER`` of -99 that the index gives an image
taken before the tour. Until the indexes are corrected, three declared values that are
real are kept: an ``EXPOSURE_DURATION`` of 0, an ``ON_CHIP_MOSAIC_FLAG`` of ``N``, and a
``PRODUCT_VERSION_ID`` of 0, which the product labels also give; any other value
declared in one of those columns is left out. Also left out are the
``SMEAR_MAGNITUDE`` of 0.1 that the index gives every image of the Jupiter tour, the
``SOLAR_DISTANCE`` of 777909100 km that the archive documents as not implemented, and
the values that it documents as not available: the ``HORIZONTAL_PIXEL_SCALE`` and
``VERTICAL_PIXEL_SCALE`` of 0 and the time from closest approach of ``-73T69:12:71Z``
that some rows give. The indexes write -999 for the other values that the labels write
as 0 where they are not available.

The current SSI indexes write each image's three compression ratios in rotated columns:
the value in ``ENCODING_COMPRESSION_RATIO`` is the image's maximum ratio, the value in
``ENCODING_MIN_COMPRESSION_RATIO`` is its ratio, and the value in
``ENCODING_MAX_COMPRESSION_RATIO`` is its minimum. For example, the product label of
``GO_0021/C21/IO/C0506406074R`` gives a ratio of 3.758, a minimum of 3.03, and a maximum
of 5.428, while the index gives 5.428, 3.758, and 3.03. Only the indexes of GO_0017
through GO_0023, and the cumulative GO_0999, give nonzero ratios. Until the indexes are
corrected, the host gives each value under the name of the ratio it holds, as the product
label does. It decides this once for each index, from the order in which its PDS3 label
lists the three columns: every SSI index label lists them in the same order, in every
version of the metadata tree, however recently it was made. An index whose label lists
the columns in another order, as a corrected one that renames them does, is taken to be
corrected, and its values keep their column names.

The same indexes give every image of the Jupiter tour the data set version 1.0, as in
``GO-J/JSA-SSI-2-REDR-V1.0``, while every product label gives version 1.1. The host
gives ``DATA_SET_ID`` as the product label does. The indexes name the time from closest
approach to the satellite ``NTV_SAT_TIME_FROM_CLOSEST_APR``, and the host gives it
under the product labels' name, ``NTV_SAT_TIME_FROM_CLOSEST_APRH``.

The subfields are those of an image read from its label (see
:doc:`user_guide_metadata`), except ``data``, ``label_dict``, and ``header_dict``. The
telemetry format comes from the supplement, so a row passed in ``row_dicts`` without the
supplement's columns gets the full field of view.

Every row gives a Snapshot, including a row whose exposure is zero and a row whose
pointing the C kernels lack. The pointing is not checked when the Snapshot is built. If
the C kernels lack the pointing of an image, its Snapshot's frame raises
:exc:`oops.OopsSpiceError` whenever it is evaluated, with the same message that
:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` gives (see :doc:`user_guide_errors`).

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
only when the image file is local or already in the local cache. Reading an index leaves
the file system as it is, even when the volumes are not there.

Reading an index through oops.Host
==================================

:meth:`oops.host.Host.from_index` recognizes an SSI index from its label and passes it,
with every other argument, to :meth:`~host_galileo.oops.ssi.GalileoSSI.from_index`,
which reads and pairs the two tables itself. The two give the same Snapshots. The label
of an SSI volume's own index does not name its instrument, but its table has the
``MISSION_NAME``, ``INSTRUMENT_ID``, and ``OBSTRUCTION_ID`` columns of an SSI volume's
index, so :meth:`oops.host.Host.from_index` recognizes it from its label too and passes
it to the host, which raises :exc:`~oops.host.HostIndexError` as described above before
any row is read. An index whose label gives neither is recognized from its first row.
