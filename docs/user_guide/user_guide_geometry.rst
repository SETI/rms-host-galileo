=============================
Time, Pointing, and Targets
=============================

Time
====

``IMAGE_TIME`` is the middle of the exposure, in UTC at the spacecraft. The exposure
begins half of ``EXPOSURE_DURATION`` earlier and ends half of it later, and the
Snapshot's default cadence covers exactly that interval.

A few images have no ``IMAGE_TIME`` in either label: 77 calibration frames in GO_0002 and
GO_0003. For these the host computes the start of the exposure from
``SPACECRAFT_CLOCK_START_COUNT``, which every label gives, and sets the
``time_from_sclk`` subfield to True:

- The SPICE clock kernel converts the count to the start of the frame in which the image
  was read out.
- The shutter opens a fixed delay after the start of that frame. The delay depends on
  ``FRAME_DURATION``: 0.225 s for frames of 2.333 s, 1.149 s for 8.667 s, 1.150 s for
  15.167 s, 2.817 s for 30.333 s, and 6.483 s for 60.667 s.
- An ``EXTENDED`` exposure taken under the phase-1 flight software opened its shutter one
  frame earlier, so the host subtracts one frame duration. The phase-2 software took
  over early in the Jupiter tour; its first image in the archive is
  ``GO_0017/J0/OPNAV/C0346405900R``.

Each delay is the median, over the 14,585 images that give both values, of the time
from the start of the frame to the start of the exposure that ``IMAGE_TIME`` implies. A
start computed this way agrees with ``IMAGE_TIME`` to within 1/3 s for about 90% of
those images.

An image cannot be read if it has neither ``IMAGE_TIME`` nor both the clock count and
the frame duration (see :doc:`user_guide_errors`).

Pointing
========

The camera's pointing comes from the scan platform C kernels, through the SPICE frame
``GLL_SCAN_PLATFORM``. The kernels sample the pointing at discrete clock ticks, and the
host accepts a sample within 40 ticks of the time it needs.

When :meth:`~host_galileo.oops.ssi.GalileoSSI.from_file` builds the frame from the C
kernels, it checks at once that they hold the pointing at the middle of the exposure,
after any ``cadence`` override. If they do not, it raises :exc:`oops.OopsSpiceError`
(see :doc:`user_guide_errors`). An index row's Snapshot is built without that check, and
its frame raises the same error only when it is evaluated (see :doc:`user_guide_index`).

A caller that has better pointing can pass it with the ``frame`` override, for example
as a :class:`~oops.frame.Cmatrix`. The host then uses that frame, does not check the C
kernels, and raises no error about missing pointing. The C kernels are still loaded,
because they belong to the fixed set of kernels that the host always furnishes.

Targets and light sources
=========================

The target body comes from ``TARGET_NAME``. A few names that the labels use are not
bodies, and the host repairs them:

.. list-table::
   :header-rows: 1

   * - ``TARGET_NAME``
     - Target
   * - ``BLACK_SKY``, ``BLACK SKY``, ``NON_SCIENCE``, ``NON SCIENCE``, ``STAR``
     - ``NONE``
   * - ``HD151288``
     - ``HD 151288``
   * - ``J RINGS``
     - ``JUPITER_RING_PLANE``

By default the light source is the Sun, except in an image of a star. When the target is
``HD 151288``, ``PLEIADES``, ``SIRIUS``, or a star name in the form ``BET_ARI``, that
star becomes the light source and the target becomes ``NONE``. ``HD 151288`` is not in
the default star catalog of OOPS, so the host adds it.

The ``label_target`` subfield keeps ``TARGET_NAME`` as the label spells it.
