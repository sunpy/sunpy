# -*- coding: utf-8 -*-

from sunpy.net.vso import VSOClient
# Import and register other sources
from sunpy.net.jsoc.jsoc import JSOCClient

# Add the JSOC and VSO Clients explicitly as they do not inherit from
# GenericClient
from .client import CLIENTS
# Import and register LC sources
from .sources.eve import EVEClient
from .sources.goes import XRSClient
from .sources.lyra import LYRAClient
from .sources.noaa import SRSClient, NOAAIndicesClient, NOAAPredictClient
from .sources.norh import NoRHClient
from .sources.rhessi import RHESSIClient

CLIENTS[VSOClient] = VSOClient._can_handle_query
CLIENTS[JSOCClient] = JSOCClient._can_handle_query
