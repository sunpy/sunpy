
# Import and register the clients but we do not want them in the namespace, we import them as _
# check if user has installed the net extras
from sunpy.util.sysinfo import warn_missing_deps
warn_missing_deps('net')

from sunpy.net import base_client as _
from sunpy.net import cdaweb as _
from sunpy.net import dataretriever as _
from sunpy.net import hek as _
from sunpy.net import helio as _
from sunpy.net import jsoc as _
from sunpy.net import vso as _
from sunpy.net.fido_factory import Fido
from sunpy.net.scraper import Scraper

__all__ = ["Fido", "Scraper"]
