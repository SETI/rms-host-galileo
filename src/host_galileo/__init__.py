##########################################################################################
# host_galileo/__init__.py
##########################################################################################
"""The Galileo mission hosts for the OOPS observation geometry toolkit.

The package exports :class:`~host_galileo.oops.ssi.GalileoSSI`. Importing
``host_galileo`` does not import oops; the first access to the class imports
:mod:`host_galileo.oops`, which registers the host with :class:`~oops.host.Host`.
"""

import importlib

__all__ = ['GalileoSSI']

# The oops hosts are imported on first access, so that importing host_galileo does not
# import oops. Importing host_galileo.oops registers every host with oops.Host.
_LAZY_NAMES = {
    'GalileoSSI': 'host_galileo.oops',
}


def __getattr__(name):
    """Import an oops host class on first access (PEP 562)."""

    if name in _LAZY_NAMES:
        return getattr(importlib.import_module(_LAZY_NAMES[name]), name)
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

##########################################################################################
