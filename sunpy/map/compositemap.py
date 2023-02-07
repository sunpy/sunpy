import sunpy.map
from sunpy.util.decorators import deprecated

__all__ = ['CompositeMap']


@deprecated("5.0", alternative="sunpy.map.CompositeMap")
class CompositeMap(sunpy.map.CompositeMap):
    pass
