"""
``sunpy-soar``
==============

A sunpy plugin for accessing data in the Solar Orbiter Archive (SOAR).

* Homepage: https://sunpy.org
* Documentation: https://docs.sunpy.org/projects/soar/
* Source Code: https://github.com/sunpy/sunpy-soar
"""

from sunpy.net.soar.attrs import SOOP, Product
from sunpy.net.soar.client import SOARClient  # Import here to register the client with sunpy

from .version import version as __version__

__all__ = ["__version__", "SOARClient", "Product", "SOOP"]
