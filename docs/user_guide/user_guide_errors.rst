======
Errors
======

The host reports each kind of failure with one exception class. The classes come from
OOPS, so code that reads images from several missions can catch them the same way.

Reading an image file
=====================

:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` raises:

:exc:`oops.OopsSpiceError`
    The C kernels lack the pointing of the scan platform at the middle of the exposure,
    after any ``cadence`` override and any ``timeshift`` of the cadence. The message
    reads "no Galileo SSI pointing at <time> for <image file>: <code>", ending with the
    SPICE error code, such as ``SPICE(CKINSUFFDATA)``. It is never raised when the
    caller supplies the pointing with the ``frame`` override. The same class is raised
    when the time must come from the clock count and the clock kernel cannot convert it;
    that message begins "no time for Galileo clock count".

:exc:`~oops.host.HostFormatError`
    No label read gives the exposure duration, the telemetry format, or a time
    (``IMAGE_TIME``, or the clock count and the frame duration); the telemetry format is
    not one of the camera's; or the time must come from the clock count and the frame
    duration is not one of the camera's. The message gives the item and the file. It
    is also raised for a normal read of a parsed label or header that has no file path,
    with a message that says no file path is known.

:exc:`~oops.host.PdsHostError`
    None of the fast, strict, and loose methods of the parser can read the PDS3 label,
    or the label has no ``^IMAGE`` pointer.

:exc:`~oops.host.VicarHostError`
    The VICAR header or image cannot be read, including an image file that is shorter
    than its header says.

:exc:`~oops.host.WrongHostError`
    The file, or the parsed label or image passed in, is not of a Galileo SSI image.

:exc:`FileNotFoundError`
    The image file that a PDS3 label points to does not exist, in a read that needs it.

:exc:`~oops.OopsTypeError`
    A keyword argument is not one that the method accepts.

:exc:`~oops.OopsValueError`
    ``select`` is not None, or ``pds3_options`` holds an option that the parser does not
    accept.

:exc:`OSError`
    A file cannot be read.

Reading a PDS3 index
====================

:meth:`~host_galileo.oops.ssi.GalileoSSI.from_index` raises these errors when it is
called:

:exc:`~oops.host.HostIndexError`
    The index is a volume's own index, in a directory named ``INDEX`` in upper or lower
    case; ``supplement`` is None and ``row_dicts`` is not given; the index's file name
    lacks the ``_index`` from which to find the supplement; a label does not describe a
    table that can be read; or a table lacks ``VOLUME_ID`` or
    ``FILE_SPECIFICATION_NAME``.

:exc:`~oops.host.PdsHostError`
    The label of the index or of the supplement cannot be parsed.

:exc:`FileNotFoundError`
    The label of the supplemental index does not exist.

:exc:`~oops.OopsTypeError`
    A keyword argument is not one that the method accepts.

:exc:`~oops.OopsValueError`
    ``select`` is not None, ``filter`` holds a column that neither table has, or
    ``pds3_options`` holds an option that the parser does not accept.

:exc:`~oops.host.WrongHostError`
    The file is not a Galileo SSI index.

:exc:`OSError`
    The label of the index or of a named supplement does not exist or cannot be read.

A table file that does not exist or cannot be read raises :exc:`OSError` when the
iteration starts. These errors are raised when the iteration reaches the row, and the
iteration ends there:

:exc:`~oops.host.HostIndexError`
    The index row and the supplement row in the same position give different images; the
    supplement ends before the index; or the table file lacks a row that its label
    promises, or a row does not match the label. The message gives the row, counted from
    0, and both tables. A supplement with rows left over raises the same error after the
    last row of the index.

:exc:`~oops.host.HostFormatError`
    A row lacks the exposure duration or a time (``IMAGE_TIME``, or the clock count and
    the frame duration), or its time must come from the clock count and its frame
    duration is not one of the camera's.

:exc:`oops.OopsSpiceError`
    A row's time must come from the clock count, and the clock kernel cannot convert it.

:exc:`~oops.OopsTypeError`
    A range in ``filter`` cannot be compared with the row's value.

Building a row's Snapshot does not check the pointing. When the C kernels lack the
pointing of an image, its Snapshot's frame raises :exc:`oops.OopsSpiceError`, with the
message described above, each time the frame is evaluated.
