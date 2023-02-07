import sunpy.data.data_manager
from sunpy.util.decorators import deprecated

__all__ = ['Cache']


@deprecated("5.0", alternative="sunpy.data.data_manager.Cache")
class Cache(sunpy.data.data_manager.Cache):
    pass
