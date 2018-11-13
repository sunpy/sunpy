__all__ = ["Fido"]


# Import and register the clients but we do not want them in the namespace, we import them as _
# Add the JSOC and VSO Clients explicitly as they do not inherit from {Base,Generic}Client
from sunpy.net import base_client as _base_client
from sunpy.net.jsoc.jsoc import JSOCClient as _
_base_client.BaseClient._registry[_] = _._can_handle_query
from sunpy.net.vso import VSOClient as _
_base_client.BaseClient._registry[_] = _._can_handle_query
from sunpy.net import dataretriever as _

from sunpy.net.fido_factory import Fido
