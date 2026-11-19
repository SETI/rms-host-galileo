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
    raised when the caller supplies the pointing with the ``frame`` override.

:exc:`~oops.host.HostFormatError`
    No label read gives the exposure duration, the telemetry format, or a time
    (``IMAGE_TIME``, or the clock count and the frame duration); the telemetry format is
    not one of the camera's; or the time must come from the clock count and the frame
    duration is not one of the camera's. The message gives the item and the file.

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

:exc:`~oops.OopsValueError`
    ``pds3_options`` holds an option that the parser does not accept.

:exc:`OSError`
    A file cannot be read.

Reading a PDS3 index
====================

:meth:`~host_galileo.oops.ssi.GalileoSSI.from_index` raises:

:exc:`~oops.host.HostIndexError`
    ``supplement`` is None and ``row_dicts`` is not given; an image appears in only one
    of the index and the supplement; the two tables list the images in different orders
    and one of them lists an image twice; or the index's file name holds no ``_index``
    from which to find the supplement.

:exc:`FileNotFoundError`
    The supplemental index does not exist.

:exc:`~oops.host.HostFormatError`
    A row gives no exposure duration or no time (``IMAGE_TIME``, or the clock count and
    the frame duration), or its time must come from the clock count and its frame
    duration is not one of the camera's. This error is raised when the iteration reaches
    that row.

:exc:`~oops.OopsTypeError`
    A keyword argument is not one that the method accepts.

:exc:`~oops.OopsValueError`
    ``select`` is not None.

:exc:`~oops.host.WrongHostError`
    The file is not a Galileo SSI index.

:exc:`OSError`
    The index or a named supplement does not exist or cannot be read.

The method does not check the pointing and never raises :exc:`oops.OopsSpiceError`.
