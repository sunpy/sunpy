import sunpy.time
from sunpy.util.decorators import deprecated

__all__ = ['TimeRange']


@deprecated("5.0", alternative="sunpy.time.TimeRange")
class TimeRange(sunpy.time.TimeRange):
    pass
