"""
The `sunpy.net.dataretriever` submodule is a framework for downloading data
from "simple" web sources such as HTTP or FTP servers. Although it could be
used for more complex services as well. Following the example of `sunpy.map`
and `sunpy.timeseries` this module provides a base class
`~sunpy.net.dataretriever.GenericClient` from which specific services can
subclass. All these subclasses are then registered with the `sunpy.net.Fido` factory
class, so do not need to be called individually.
"""

from .client import GenericClient, QueryResponse
from .sources.eve import EVEClient
from .sources.fermi_gbm import GBMClient
from .sources.goes import SUVIClient, XRSClient
from .sources.gong import GONGClient
from .sources.lyra import LYRAClient
from .sources.noaa import NOAAIndicesClient, NOAAPredictClient, SRSClient
from .sources.norh import NoRHClient
from .sources.rhessi import RHESSIClient

__all__ = ['QueryResponse', 'GenericClient', 'EVEClient', 'GBMClient', 'XRSClient',
           'GONGClient', 'SUVIClient', 'LYRAClient', 'NOAAIndicesClient',
           'NOAAPredictClient', 'NoRHClient', 'RHESSIClient', 'SRSClient']
