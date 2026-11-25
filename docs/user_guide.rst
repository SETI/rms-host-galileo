User Guide
==========

``rms-host-galileo`` turns Galileo Solid State Imaging (SSI) images into OOPS
:class:`~oops.observation.Snapshot` objects. Importing ``host_galileo.oops`` registers the
host with :class:`oops.Host <oops.host.Host>`, so that
:meth:`oops.Host.from_file <oops.host.Host.from_file>` and
:meth:`oops.Host.from_index <oops.host.Host.from_index>` recognize SSI files. The host class,
:class:`~host_galileo.oops.ssi.GalileoSSI`, can also be used directly:

.. code-block:: python

    from host_galileo import GalileoSSI

    obs = GalileoSSI.from_file('GO_0017/G1/GANYMEDE/C0349632100R.LBL')
    for obs in GalileoSSI.from_index('GO_0017_index.lbl'):
        ...


Setup
-----

Geometry needs the Galileo SPICE kernels. The host reads a fixed set of them through the
SPICE file cache of ``spicedb``, located by the ``SPICE_PATH`` environment variable (or
``$OOPS_RESOURCES/SPICE``), and the SSI instrument kernel through the kernel database
named by ``SPICE_SQLITE_DB_NAME``.

Importing the host does not load any kernels. The first call to
:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` or
:meth:`~host_galileo.oops.ssi.GalileoSSI.from_index` calls
:meth:`~host_galileo.oops.ssi.GalileoSSI.initialize` with its defaults, which defines the
solar system bodies and loads the kernels. Call it yourself first, with keyword
arguments, to choose other options: ``planets`` limits the planets that are defined,
``asof`` uses only the kernels that existed before a date, and ``irregulars=False``
leaves out the irregular moons. Later calls change nothing.


Reading an image
----------------

:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` takes either file of an SSI product:
the detached PDS3 label (a name ending in ``.lbl``, in any case) or the VICAR image file.

- Given the image file, only that file is read.
- Given the label, the host also reads the image file that the label's ``^IMAGE`` pointer
  names, in the label's directory and in whatever case it has there.
- With ``astrometry=True``, the image file of a label is not read at all, and no image
  array is read from any file. Given the image file, its header is still read.

The image array always comes from the VICAR file, laid out as its own header describes.
The record number in the label's ``^IMAGE`` pointer is not used, because some labels give
it wrongly: 52 labels of GO_0019 put the image one record early, and the label of GO_0018
E4/EUROPA/C0374685452R describes the wrong layout. One read of the image file gives both
its header and its array.

The cutout window
~~~~~~~~~~~~~~~~~

Many SSI images hold data only inside a cutout window, given by ``CUT_OUT_WINDOW`` as the
first line, first sample, number of lines, and number of samples, counted from one. By
default the array is trimmed to that window and the field of view is the matching slice
of the camera's. With ``full_fov=True``, the full 800x800 array and field of view are
returned instead. Either way, the ``window`` subfield gives the window as zero-based
line and sample limits, ``[l0, s0, l1, s1]``, or ``[0, 0, 800, 800]`` for an image with
no window.

Time
~~~~

``IMAGE_TIME`` is the middle of the exposure. Some calibration frames of GO_0002 and
GO_0003 have no time: the label writes ``IMAGE_TIME = UNK``, and the header writes either
-32768 as the year (48 frames, all in GO_0002) or 0 for all six time fields (29 frames,
22 in GO_0002 and 7 in GO_0003). Their exposure starts at the time
:meth:`~host_galileo.oops.ssi.GalileoSSI.tstart_from_sclk` derives from the spacecraft
clock count, and their ``time_from_sclk`` subfield is True. An exposure of zero is kept
as zero.

Target and light source
~~~~~~~~~~~~~~~~~~~~~~~

The target is the label's ``TARGET_NAME``, with a few names replaced by the body that
OOPS knows: ``BLACK_SKY``, ``BLACK SKY``, ``NON_SCIENCE``, ``NON SCIENCE``, and ``STAR``
become ``NONE``, ``J RINGS`` becomes ``JUPITER_RING_PLANE``, and ``HD151288`` becomes
``HD 151288``. An image of a star has no target body; the star is its light source.
Every other image is lit by the Sun. The ``label_target`` subfield keeps the name as the
label gives it.

Pointing
~~~~~~~~

The camera frame is ``GLL_SCAN_PLATFORM``, whose orientation comes from the Galileo C
kernels. The host does not check the kernels when it builds an observation. If they lack
the pointing for an image, the SPICE error is raised when the frame is first used, and
it names the SPICE error code. A ``frame`` override, such as a fixed
:class:`~oops.frame.Cmatrix`, is used as given.

Choosing from the result
~~~~~~~~~~~~~~~~~~~~~~~~

An SSI image holds one observation. The ``select`` option chooses from it as from a
tuple holding that one observation: ``0`` or ``-1`` gives the observation, and a slice
gives a tuple. Without ``select``, the result is the observation itself.


Metadata
--------

Each observation has three dictionaries of metadata:

- ``dict``: the merged metadata of every label that was read.
- ``label_dict``: the PDS3 label, exactly as the PDS3 parser gives it, if the label was
  read. It includes the keys the parser adds beside each keyword, such as
  ``IMAGE_TIME_fmt``.
- ``header_dict``: the VICAR header, exactly as the VICAR reader gives it, if the image
  file was read. A parameter that occurs more than once is keyed by its name and its
  occurrence, such as ``('TASK', 1)``. Its ``HOST`` and ``BHOST`` are those of the
  machine that reads the file, because the VICAR reader replaces the values in the file.

How ``dict`` is built
~~~~~~~~~~~~~~~~~~~~~

The PDS3 label gives its values under their own names. A header value that the host
converts to a PDS3 item appears only under the PDS3 name, with the type the label gives
it:

- the keywords that name the same value, such as ``EXP`` for ``EXPOSURE_DURATION``,
  ``PA`` for ``OBSERVATION_ID``, ``TARGET`` for ``TARGET_NAME``, ``FILTER`` for
  ``FILTER_NUMBER``, and ``TLMFMT`` for ``TELEMETRY_FORMAT_ID``;
- ``RIM`` and ``MOD91`` as ``SPACECRAFT_CLOCK_START_COUNT``, written ``RRRRRRRR.MM``;
- the ``SCET`` and ``ERT`` fields as ``IMAGE_TIME`` and ``EARTH_RECEIVED_TIME``, each a
  :class:`datetime.datetime`, as the label's times are;
- ``RATE`` as ``FRAME_DURATION`` and ``GAIN`` as ``GAIN_MODE_ID``;
- ``FILTER_NAME`` from the filter number, and ``EXPOSURE_TYPE`` from the last of the
  ``MOFIBE`` flags (or the ``FIBE`` flags of the frames taken before the Jupiter tour).

``MOD10`` and ``MOD8``, the finer fields of the clock count, keep their own names, since
no PDS3 item holds them. Every other header value keeps its own name too.

An item is left out of ``dict`` when no label read gives it, or when its value means
"no value": a PDS3 null text (``UNK``, ``N/A``, ``NULL``, ``UNKNOWN``,
``NOT APPLICABLE``, ``---``, or an empty string, in any case), or a header value of
-999. A header that writes -32768 or 0 as the year gives no ``IMAGE_TIME``. Every other value
is given as the file writes it.

A label's ``ORBIT_NUMBER`` outside 0 to 35, the orbits of the Jupiter tour, is left out.
The labels of GO_0003 to GO_0009 and GO_0016 give the frames taken before the tour an
``ORBIT_NUMBER`` made of two characters of the observation ID read as a number, such as
12613 for ``E1``. The host derives no ``ORBIT_NUMBER`` from a header, so an image file
read alone has none.

Which label wins
~~~~~~~~~~~~~~~~

Where both labels give an item and the values differ, each item takes its value from one
fixed label for every image. The choice rests on a comparison of the two labels of all
14,672 SSI products, and on checks against SPICE.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Winner
     - Items
   * - VICAR header
     - ``EXPOSURE_DURATION``, which the label rounds to two decimal places; the pointing
       and viewing geometry that the label rounds or, in some frames, miscomputes:
       ``CONE_ANGLE``, ``DECLINATION``, ``EMISSION_ANGLE``, ``ENTROPY``,
       ``INCIDENCE_ANGLE``, ``INTERCEPT_POINT_LATITUDE``, ``INTERCEPT_POINT_LINE``,
       ``INTERCEPT_POINT_LINE_SAMPLE``, ``INTERCEPT_POINT_LONGITUDE``,
       ``LOCAL_HOUR_ANGLE``, ``NORTH_AZIMUTH``, ``PHASE_ANGLE``, ``RIGHT_ASCENSION``,
       ``SUB_SOLAR_AZIMUTH``, ``SUB_SOLAR_LATITUDE``, ``SUB_SOLAR_LONGITUDE``,
       ``SUB_SPACECRAFT_AZIMUTH``, ``SUB_SPACECRAFT_LATITUDE``,
       ``SUB_SPACECRAFT_LINE``, ``SUB_SPACECRAFT_LINE_SAMPLE``,
       ``SUB_SPACECRAFT_LONGITUDE``, ``TARGET_CENTER_DISTANCE``, and ``TWIST_ANGLE``;
       and the compression ratios, ``ENCODING_COMPRESSION_RATIO``,
       ``ENCODING_MAX_COMPRESSION_RATIO``, and ``ENCODING_MIN_COMPRESSION_RATIO``.
   * - PDS3 label
     - Every other item, including ``CENTRAL_BODY_DISTANCE``, ``COMPRESSION_TYPE``, and
       ``TARGET_NAME``.

If the winning label lacks the item, the other label's value is used.


Subfields
---------

Besides ``dict``, ``label_dict``, and ``header_dict``, each observation from
:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` has these subfields:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Subfield
     - Value
   * - ``data``
     - The image array, unless ``astrometry`` is True.
   * - ``host``, ``instrument``
     - ``"Galileo"`` and ``"SSI"``.
   * - ``filter``, ``label_target``, ``observation_id``, ``spacecraft_clock_start_count``
     - ``FILTER_NAME``, ``TARGET_NAME``, ``OBSERVATION_ID``, and
       ``SPACECRAFT_CLOCK_START_COUNT`` from ``dict``, each only if ``dict`` has it.
   * - ``time_from_sclk``
     - True if the time comes from the spacecraft clock count.
   * - ``window``, ``full_fov``
     - The cutout window, as described above, and the ``full_fov`` option.
   * - ``spice_kernels``
     - The names of the SPICE kernels that apply to the observation. With a ``frame``
       override, no C kernel is listed.
   * - ``spice_frame_name``, ``spice_frame_id``, ``spice_to_frame``
     - ``"GLL_SCAN_PLATFORM"``, -77001, and the identity rotation from that SPICE frame
       to the camera frame.
   * - ``ck_object_id``, ``ck_sclk_id``, ``spice_frame_is_frozen``
     - -77001 and -77, the NAIF IDs of the object whose C kernels give the pointing and
       of the spacecraft clock they use, and False, since the frame follows the C
       kernels.
   * - ``filepath``, ``basename``, ``image_url``
     - The image file, also when the label was given: its path, its name, and its
       absolute path or URL.
   * - ``image_abspath``
     - The path of the image file on the local disk, if the host read the image file.
   * - ``label_url``, ``label_abspath``
     - The label's absolute path or URL and its local path, if the label was read.


Reading an index
----------------

:meth:`~host_galileo.oops.ssi.GalileoSSI.from_index` reads a metadata index of the PDS
Ring-Moon Systems Node, such as ``GO_0017_index.lbl``, together with its supplemental
index, and returns an iterator that yields one observation for each row. Each
observation is built when the iteration reaches its row, and the tables are read 100
rows at a time. A caller that needs ``len()`` or indexing wraps the call in ``list()``.

- The supplemental index is required, because only it gives ``TELEMETRY_FORMAT_ID`` and
  ``CUT_OUT_WINDOW``. By default it is the file beside the index whose name has
  ``_supplemental_index`` in place of ``_index``. A ``supplement`` path in any directory
  is read instead. ``supplement=None`` raises an error, and so does a supplement that
  does not exist.
- Every SSI index lists its products in the same order as its supplement, so each index
  row is paired with the supplement row in the same position, and both must name the
  same file. No other column is compared, because the two tables write
  ``SPACECRAFT_CLOCK_START_COUNT`` in different forms (``03496321.00`` in the index,
  ``03496321:00:0:0`` in the supplement). Where both give a column, the index value is
  kept.
- ``dict`` is the merged row, with every value as the table reader gives it, declared
  null values included, along with the reader's ``_mask`` entry for each column.
- Every row gives an observation, including the rows whose exposure is zero.
- A ``CUT_OUT_WINDOW`` of ``-1, -1, -1, -1`` means that the image has no cutout window.
- An ``IMAGE_TIME`` of ``UNK`` means that the time comes from the spacecraft clock
  count.
- ``filepath`` is the row's ``VOLUME_ID`` and ``FILE_SPECIFICATION_NAME`` joined by a
  slash, such as ``GO_0017/G1/GANYMEDE/C0349632100R.LBL``, and ``basename`` is its last
  part. To reach the file, join ``filepath`` to your own copy of the volumes. A row has
  no other file subfield, no ``data``, and no ``label_dict`` or ``header_dict``.

The ``filter`` option keeps only the rows whose columns have the given values, and
``row_dicts`` takes rows that have already been read and merged, in place of the tables.


The start time from the clock count
-----------------------------------

:meth:`~host_galileo.oops.ssi.GalileoSSI.tstart_from_sclk` takes a label as a dictionary,
or a merged index row, and returns the start of the exposure in seconds TDB. It needs
``SPACECRAFT_CLOCK_START_COUNT``, in the label form ``RRRRRRRR.MM`` or the SPICE form
``RRRRRRRR:MM:S:T``, and ``FRAME_DURATION``; ``EXPOSURE_TYPE`` is taken as ``NORMAL`` if
it is absent. It loads only the leap-second and spacecraft clock kernels, on its first
call, and needs no call to :meth:`~host_galileo.oops.ssi.GalileoSSI.initialize`.

The clock count names the frame in which the image was read out. The shutter opens a
fixed delay after the start of that frame, which depends on the frame duration: 0.225 s
for 2.333 s frames, 1.149 s for 8.667 s, 1.150 s for 15.167 s, 2.817 s for 30.333 s, and
6.483 s for 60.667 s. An ``EXTENDED`` exposure taken before the Jupiter tour's flight
software was updated, which is any frame before GO_0017 J0/OPNAV/C0346405900R (clock
count 03464059.00), opens one frame earlier. The result agrees with the start that
``IMAGE_TIME`` implies to within about 1/3 s for most images.


Errors
------

Every failure to read or interpret a label or an image file raises an error from
:mod:`oops.host` that names the file, with one exception: when ``oops.Host.from_file``
identifies the file as a VICAR image, an image file that is shorter than its header says
raises a :class:`ValueError` instead of ``VicarHostError``. The same file read through
:meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` raises ``VicarHostError``.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Error
     - Raised when
   * - ``PdsHostError``
     - The PDS3 label cannot be parsed, including a value the parser cannot convert, or
       it has no ``^IMAGE`` pointer.
   * - ``VicarHostError``
     - The image file is not a valid VICAR file, is shorter than its header says, or
       holds more than one band.
   * - ``HostFormatError``
     - The array does not have the dimensions that the label gives it, or the metadata
       cannot be used: for example, no label gives the exposure, the telemetry format,
       or a time, the telemetry format is not one of SSI's, or the frame duration is
       not one of SSI's when the time must come from the clock count.
   * - ``WrongHostError``
     - The file is not that of a Galileo SSI image, or the index is not an SSI index.
   * - ``HostIndexError``
     - ``supplement`` is None, the supplement's name cannot be derived from the index's,
       or the rows of the two tables do not name the same files.

A missing file raises :class:`FileNotFoundError`, and an unknown keyword argument raises
:class:`TypeError`. The rows of an index are used as the table gives them, with no
check.
