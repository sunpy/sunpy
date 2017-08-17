"""
The `sunpy.net.dataretriever` submodule is a framework for downloading data
from "simple" web sources such as HTTP or FTP servers. Although it could be
used for more complex services as well. Following the example of `sunpy.map`
and `sunpy.timeseries` this module provides a base class
`~sunpy.net.dataretriever.GenericClient` from which specific services can
subclass. All these subclasses are then registered with the `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` factory
class, so do not need to be called individually.
"""

from .client import QueryResponseBlock, QueryResponse, GenericClient

from . import clients
