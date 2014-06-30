from .factbase import UnifiedDownloader
from .eve import EVEClient
from .lyra import LYRAClient
from .norh import NoRHClient
from ..vso import VSOClient
UnifiedDownloader.register(EVEClient,EVEClient._can_handle_query)
UnifiedDownloader.register(LYRAClient,LYRAClient._can_handle_query)
UnifiedDownloader.register(NoRHClient,NoRHClient._can_handle_query)
UnifiedDownloader.register(VSOClient,VSOClient._can_handle_query)

