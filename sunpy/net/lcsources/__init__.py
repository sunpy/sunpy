from .factbase import UnifiedDownloader
from .eve import EVEClient
UnifiedDownloader.register(EVEClient,EVEClient._can_handle_query)

