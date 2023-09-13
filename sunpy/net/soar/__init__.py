# Import here to register the client with sunpy
from sunpy.net.soar.attrs import SOOP, Product
from sunpy.net.soar.client import SOARClient

from .version import version as __version__

__all__ = ["__version__", "SOARClient", "Product", "SOOP"]
