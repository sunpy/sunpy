import sunpy.map
from sunpy.util.decorators import deprecated

__all__ = ['MapSequence']


@deprecated("5.0", alternative="sunpy.map.MapSequence")
class MapSequence(sunpy.map.MapSequence):
    pass
