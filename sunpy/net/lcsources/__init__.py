from .factbase import UnifiedDownloader
from .eve import EVEClient
from .lyra import LYRAClient
from .norh import NoRHClient
UnifiedDownloader.register(EVEClient,EVEClient._can_handle_query)
UnifiedDownloader.register(LYRAClient,LYRAClient._can_handle_query)
UnifiedDownloader.register(NoRHClient,NoRHClient._can_handle_query)

