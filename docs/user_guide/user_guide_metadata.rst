=======================
Metadata and Subfields
=======================

The two labels
==============

An SSI product describes itself twice. The *PDS3 label* is the detached ``.LBL`` file
beside the image. The *VICAR header* is the label at the start of the image file itself.
The two were written at different times and do not always agree, and each holds values
that the other lacks.

The host takes its metadata from whichever of the two it read (see
:doc:`user_guide_files`): from the VICAR header alone when given the image file, from
the PDS3 label alone in an astrometry read of the label, and from both when given the
label otherwise.

Merging the two labels
======================

The VICAR header uses its own keywords and codes. Where a header value has a PDS3
equivalent, the host converts it to the PDS3 name and to the form that the PDS3 label
uses. For example:

- ``EXP`` becomes ``EXPOSURE_DURATION``, ``PA`` becomes ``OBSERVATION_ID``, ``TARGET``
  becomes ``TARGET_NAME``, and ``TLMFMT`` becomes ``TELEMETRY_FORMAT_ID``.
- ``RIM`` and ``MOD91`` together become ``SPACECRAFT_CLOCK_START_COUNT``.
- The ``SCET`` and ``ERT`` time fields become ``IMAGE_TIME`` and
  ``EARTH_RECEIVED_TIME``, each written as ISO 8601 text such as
  ``1996-06-26T08:46:04.890``.
- The codes ``RATE`` and ``GAIN`` become ``FRAME_DURATION`` and ``GAIN_MODE_ID``, and
  the filter number gives ``FILTER_NAME``.
- The flags ``FIBE`` or ``MOFIBE`` give ``EXPOSURE_TYPE``, and an observation ID that
  begins with a Jupiter tour orbit gives ``ORBIT_NUMBER``.

A converted value appears in the merged metadata only under its PDS3 name, and the
original keyword stays in ``header_dict``. A header keyword that has no PDS3 equivalent
keeps its own name. ``MOD10`` and ``MOD8``, the finer fields of the clock count, are
examples, because the PDS3 count holds only the ``RIM`` and ``MOD91`` fields.

Where only one label gives an item, its value is used. Where both give it and agree,
that value is used. Where both give it and disagree, the item always takes its value
from the same label, for every image. The VICAR header supplies these items:

.. list-table::
   :header-rows: 1

   * - Items
     - Why the VICAR header
   * - ``EXPOSURE_DURATION``
     - The header gives each commanded exposure, such as 4 1/6 ms, to within 1 µs. The
       label rounds it to two decimal places.
   * - ``CONE_ANGLE``, ``DECLINATION``, ``EMISSION_ANGLE``, ``ENTROPY``,
       ``INCIDENCE_ANGLE``, ``INTERCEPT_POINT_LATITUDE``, ``INTERCEPT_POINT_LINE``,
       ``INTERCEPT_POINT_LINE_SAMPLE``, ``INTERCEPT_POINT_LONGITUDE``,
       ``LOCAL_HOUR_ANGLE``, ``NORTH_AZIMUTH``, ``PHASE_ANGLE``, ``RIGHT_ASCENSION``,
       ``SUB_SOLAR_AZIMUTH``, ``SUB_SPACECRAFT_AZIMUTH``, ``SUB_SPACECRAFT_LONGITUDE``,
       and ``TWIST_ANGLE``
     - The label gives the header's value rounded to fewer digits.
   * - ``SUB_SOLAR_LATITUDE``, ``SUB_SOLAR_LONGITUDE``, ``SUB_SPACECRAFT_LATITUDE``,
       ``SUB_SPACECRAFT_LINE``, ``SUB_SPACECRAFT_LINE_SAMPLE``, and
       ``TARGET_CENTER_DISTANCE``
     - Where the two differ by more than rounding, the header agrees with the SPICE
       kernels and the label does not.
   * - ``ENCODING_COMPRESSION_RATIO``, ``ENCODING_MAX_COMPRESSION_RATIO``, and
       ``ENCODING_MIN_COMPRESSION_RATIO``
     - The label copies the ratios of an on-chip mosaic's parent image to each of its
       parts. The header gives each image's own.

The PDS3 label supplies every other item, including ``TARGET_NAME``,
``CENTRAL_BODY_DISTANCE``, and ``COMPRESSION_TYPE``. If the chosen label gives no value
for an item, the other label's value is used.

Missing values
==============

A value that no label read gives is left out: there is no ``dict`` key and no subfield
for it, and never a None or a placeholder in its place. These count as no value:

- the PDS3 null constants ``UNK``, ``N/A``, and ``NULL``, in any case, and a list of
  nothing else;
- a numeric -999 in the VICAR header, and a ``PICNO`` or ``SPICE_C_ID`` of ``?``;
- a set of VICAR time fields that makes no valid time, such as an hour of 24;
- an ``ORBIT_NUMBER`` outside the Jupiter tour's orbits 0 to 35. The labels of the
  earlier volumes give the first two characters of the observation ID read as a number,
  such as 12613 for ``E1``.

The dictionaries
================

Three subfields hold the metadata as dictionaries. All three can be written with
:func:`json.dump`: times are ISO 8601 text, sets are lists, and an item has the same
type whichever label it came from.

``dict``
    The merged metadata: every item of either label read, under the PDS3 name where the
    host converts it, and under its own name otherwise.

``label_dict``
    The PDS3 label as its parser reads it, present only when the label was read. Times
    are the label's own text, and sets are lists in the label's order. The ``END``
    keyword and the extra keys that the parser adds beside each keyword, such as
    ``IMAGE_TIME_fmt`` and ``^IMAGE_offset``, are left out.

``header_dict``
    Every parameter of the VICAR header, present only when the image file was read. A
    parameter that occurs more than once, such as the ``TASK`` of each step of the
    processing history, is keyed by its name and its occurrence counted from zero:
    ``TASK[0]``, ``TASK[1]``, and so on.

Subfields
=========

Besides ``dict``, ``label_dict``, ``header_dict``, and ``data``, every Snapshot has
these subfields. One whose item no label gives is absent.

.. list-table::
   :header-rows: 1

   * - Subfield
     - Contents
   * - ``host``, ``instrument``
     - ``"Galileo"`` and ``"SSI"``.
   * - ``filter``
     - ``FILTER_NAME``, such as ``"CLEAR"`` or ``"IR-7560"``.
   * - ``label_target``
     - ``TARGET_NAME`` as the label spells it.
   * - ``observation_id``
     - ``OBSERVATION_ID``.
   * - ``spacecraft_clock_start_count``
     - ``SPACECRAFT_CLOCK_START_COUNT``, described below.
   * - ``time_from_sclk``
     - True if no label read gives ``IMAGE_TIME``, so that the time comes from the clock
       count (see :doc:`user_guide_geometry`).
   * - ``window``, ``full_fov``
     - The cutout window and the ``full_fov`` option (see :doc:`user_guide_products`).
   * - ``spice_kernels``
     - The SPICE kernels used for the observation.
   * - ``spice_frame_name``, ``spice_frame_id``, ``spice_to_frame``
     - The SPICE frame of the camera, ``GLL_SCAN_PLATFORM`` (-77001), and the rotation
       from it to the camera frame, which is the identity.
   * - ``filepath``, ``basename``, ``abspath``, ``image_url``
     - The image file (see :doc:`user_guide_files`).

The spacecraft clock count
==========================

``spacecraft_clock_start_count`` is the count written as the PDS3 label writes it,
``RRRRRRRR.MM``, such as ``03496321.37``: the RIM count and the mod-91 count. From the
VICAR header it is formed from ``RIM`` and ``MOD91``. The header's finer fields,
``MOD10`` and ``MOD8``, are in ``dict`` and ``header_dict`` when the image file was
read. Neither label gives a stop count, so there is no
``spacecraft_clock_stop_count`` subfield.
