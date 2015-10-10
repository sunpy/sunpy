# -*- coding: utf-8 -*-

from .downloader_factory import Fido as _Fido

# Import and register LC sources
from .sources.eve import EVEClient
from .sources.lyra import LYRAClient
from .sources.goes import GOESClient
from .sources.norh import NoRHClient
from .sources.rhessi import RHESSIClient
from .sources.noaa import NOAAIndicesClient, NOAAPredictClient
from .sources.stereo import SEPTClient, HETClient, SITClient, PLASTICClient, MAGClient, LETClient
from .sources.soho import ERNEClient



_Fido.register(EVEClient, EVEClient._can_handle_query)
_Fido.register(LYRAClient, LYRAClient._can_handle_query)
_Fido.register(GOESClient, GOESClient._can_handle_query)
_Fido.register(NoRHClient, NoRHClient._can_handle_query)
_Fido.register(NOAAIndicesClient, NOAAIndicesClient._can_handle_query)
_Fido.register(NOAAPredictClient, NOAAPredictClient._can_handle_query)
_Fido.register(RHESSIClient, RHESSIClient._can_handle_query)
_Fido.register(ERNEClient, ERNEClient._can_handle_query)
_Fido.register(SEPTClient, SEPTClient._can_handle_query)
_Fido.register(HETClient, HETClient._can_handle_query)
_Fido.register(LETClient, LETClient._can_handle_query)
_Fido.register(PLASTICClient, PLASTICClient._can_handle_query)
_Fido.register(SITClient, SITClient._can_handle_query)
# _Fido.register(MAGClient, MAGClient._can_handle_query)   ** CDF data file issue ... mail sent to dev-mailing list **



# Import and register other sources
from sunpy.net.jsoc.jsoc import JSOCClient
from sunpy.net.vso.vso import VSOClient

_Fido.register(VSOClient, VSOClient._can_handle_query)
_Fido.register(JSOCClient, JSOCClient._can_handle_query)
