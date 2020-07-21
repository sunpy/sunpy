
# Import and register the clients but we do not want them in the namespace, we import them as _
from sunpy.net import base_client as _
from sunpy.net import dataretriever as _
from sunpy.net import jsoc as _
from sunpy.net import vso as _
from sunpy.net.fido_factory import Fido
from sunpy.net.hek import hek as _
from sunpy.net.helio import hec as _

__all__ = ["Fido"]
