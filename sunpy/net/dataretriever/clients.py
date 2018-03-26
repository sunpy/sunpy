# -*- coding: utf-8 -*-

# Import and register LC sources
from .sources.eve import EVEClient
from .sources.lyra import LYRAClient
from .sources.goes import XRSClient
from .sources.norh import NoRHClient
from .sources.rhessi import RHESSIClient
from .sources.noaa import NOAAIndicesClient, NOAAPredictClient, SRSClient

# Import and register other sources
from sunpy.net.jsoc.jsoc import JSOCClient
from sunpy.net.vso import VSOClient

# Add the JSOC and VSO Clients explicitly as they do not inherit from
# GenericClient
from .client import CLIENTS
CLIENTS[VSOClient] = VSOClient._can_handle_query
CLIENTS[JSOCClient] = JSOCClient._can_handle_query
