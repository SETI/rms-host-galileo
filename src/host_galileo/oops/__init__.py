##########################################################################################
# host_galileo/oops/__init__.py
##########################################################################################
"""The Galileo hosts for oops.

The host is :class:`~host_galileo.oops.ssi.GalileoSSI`, documented with its module below.
Importing this package registers it with :class:`~oops.host.Host`, so that its
``from_file()`` and ``from_index()`` recognize Galileo SSI images and indexes.
"""

# Importing each instrument module registers its host with oops.Host
from .ssi import GalileoSSI

__all__ = ['GalileoSSI']

##########################################################################################
