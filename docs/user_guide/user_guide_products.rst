==================
Supported Products
==================

Instrument
==========

The Solid State Imaging (SSI) camera of the Galileo orbiter took images of 800 by 800
pixels on a CCD mounted on the spacecraft's scan platform. The host reads the SSI images
of the PDS3 volumes GO_0002 through GO_0023, from the cruise to the end of the Jupiter
tour.

Products
========

Each SSI product is a VICAR image file, such as ``C0349632100R.IMG``, with a detached
PDS3 label beside it, such as ``C0349632100R.LBL``. The VICAR file holds a header of its
own, so a product carries two labels. :doc:`user_guide_metadata` describes how the host
combines them.

Every image becomes an OOPS :class:`~oops.observation.Snapshot`, with the image in its
``data`` subfield and the image's lines and samples as its *v* and *u* axes.

Telemetry formats and fields of view
====================================

The telemetry format, given by ``TELEMETRY_FORMAT_ID``, decides the field of view. The
full field of view is 800 by 800 pixels, with the camera's barrel distortion from the
SSI instrument kernel. The summed formats ``HIS`` and ``AI8`` combine each 2 by 2 block
of pixels into one, so their field of view is 400 by 400 pixels. The formats ``FULL``,
``HMA``, ``HCA``, ``HIM``, ``IM8``, ``IM4``, ``XCM``, ``HCJ``, ``HCM``, and ``NONE`` have
the full field of view.

Cutout windows
==============

Many images were read out only within a rectangle of the frame, given by the label's
``CUT_OUT_WINDOW``. By default the host trims the image array and the field of view to
that window. With ``full_fov=True``, both cover the full frame. The ``window`` subfield
gives the window as zero-based line and sample limits, ``[l0, s0, l1, s1]``, or the full
frame when the image has no window.

Calibration
===========

The observation's list of calibrations is empty unless the caller supplies one with the
``calibrations`` override, because the host supplies none of its own.
