======
Errors
======

The host reports each kind of failure with one exception class. The classes come from
OOPS, so code that reads images from several missions can catch them the same way.

Reading an image file
=====================

:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` raises:

:exc:`oops.OopsSpiceError`
    The C kernels hold no pointing for the scan platform at the middle of the exposure.
    The message reads "no Galileo SSI pointing at <time> for <image file>: <code>",
    ending with the SPICE error code, such as ``SPICE(CKINSUFFDATA)``. It is never
    raised when the caller supplies the pointing with the ``frame`` override. The same
    class is raised when the time must come from the clock count and the clock kernel
    cannot convert it; that message begins "no time for Galileo clock count".

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
    The index is a volume's own index, in its ``INDEX`` directory; ``supplement`` is None
    and ``row_dicts`` is not given; either table cannot be read as an index table; an
    image appears in only one of the index and the supplement; the two tables list the
    images in different orders and one of them lists an image twice; or the index's file
    name holds no ``_index`` from which to find the supplement.

:exc:`~oops.host.PdsHostError`
    The label of the index or of the supplement cannot be parsed.

:exc:`FileNotFoundError`
    The supplemental index does not exist.

:exc:`~oops.OopsTypeError`
    A keyword argument is not one that the method accepts, or a range in ``filter``
    cannot be compared with the values of its column.

:exc:`~oops.OopsValueError`
    ``select`` is not None, ``filter`` holds a column that no row has, or
    ``pds3_options`` holds an option that the parser does not accept.

:exc:`~oops.host.WrongHostError`
    The file is not a Galileo SSI index.

:exc:`OSError`
    The index or a named supplement does not exist or cannot be read.

These errors are raised when the iteration reaches the row, and the iteration ends
there:

:exc:`~oops.host.HostFormatError`
    A row gives no exposure duration or no time (``IMAGE_TIME``, or the clock count and
    the frame duration), or its time must come from the clock count and its frame
    duration is not one of the camera's.

:exc:`oops.OopsSpiceError`
    A row's time must come from the clock count, and the clock kernel cannot convert it.

Building a row's Snapshot does not check the pointing. When the C kernels hold no
pointing for an image, its Snapshot's frame raises :exc:`oops.OopsSpiceError`, with the
message described above, each time the frame is evaluated.
