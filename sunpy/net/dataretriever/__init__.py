"""
The `sunpy.net.dataretriever` submodule is a framework for downloading data
from "simple" web sources such as HTTP or FTP servers. Although it could be
used for more complex services as well. Following the example of `sunpy.map`
and `sunpy.timeseries` this module provides a base class
`~sunpy.net.dataretriever.GenericClient` from which specific services can
subclass. All these subclasses are then registered with the `sunpy.net.Fido` factory
class, so do not need to be called individually.
"""

from .client import QueryResponseBlock, QueryResponse, GenericClient
from .sources.eve import EVEClient
from .sources.lyra import LYRAClient
from .sources.goes import XRSClient, SUVIClient
from .sources.norh import NoRHClient
from .sources.rhessi import RHESSIClient
from .sources.noaa import NOAAIndicesClient, NOAAPredictClient, SRSClient
from .sources.fermi_gbm import GBMClient
from .sources.bbso import BBSOClient

__all__ = ['QueryResponseBlock', 'QueryResponse', 'GenericClient',
           'EVEClient', 'XRSClient', 'SUVIClient', 'LYRAClient', 'NOAAIndicesClient',
           'NOAAPredictClient', 'NoRHClient', 'RHESSIClient', 'SRSClient', 'GBMClient',
           'BBSOClient']
