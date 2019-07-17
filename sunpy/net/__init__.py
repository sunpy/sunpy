# -*- coding: utf-8 -*-

# Import and register the clients but we do not want them in the namespace, we import them as _
from sunpy.net import base_client as _, dataretriever as _
from sunpy.net import vso as _, jsoc as _
from sunpy.net.fido_factory import Fido

__all__ = ["Fido"]
