"""
``sunpy-soar``
==============

A sunpy FIDO plugin for accessing data in the Solar Orbiter Archive (SOAR).

* Homepage: https://sunpy.org
* Documentation: https://docs.sunpy.org/projects/soar/
* Source Code: https://github.com/sunpy/sunpy-soar
"""

# Import here to register the client with sunpy
from sunpy.net.soar.client import SOARClient

from .version import version as __version__

__all__ = ["SOARClient", "__version__"]
