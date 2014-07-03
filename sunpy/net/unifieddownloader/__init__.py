from sunpy.net.unifieddownloader.downloader_factory import UnifiedDownloader
from sunpy.net.unifieddownloader.sources.eve import EVEClient
from sunpy.net.unifieddownloader.sources.lyra import LYRAClient
from sunpy.net.unifieddownloader.sources.norh import NoRHClient
from sunpy.net.vso.vso import VSOClient
from sunpy.net.unifieddownloader.sources.rhessi import RHESSIClient
from sunpy.net.unifieddownloader.sources.noaa import NOAAIndicesClient,NOAAPredictClient
UnifiedDownloader.register(EVEClient,EVEClient._can_handle_query)
UnifiedDownloader.register(LYRAClient,LYRAClient._can_handle_query)
UnifiedDownloader.register(NoRHClient,NoRHClient._can_handle_query)
UnifiedDownloader.register(VSOClient,VSOClient._can_handle_query)
UnifiedDownloader.register(NOAAIndicesClient,NOAAIndicesClient._can_handle_query)
UnifiedDownloader.register(NOAAPredictClient,NOAAPredictClient._can_handle_query)
UnifiedDownloader.register(RHESSIClient,RHESSIClient._can_handle_query)
