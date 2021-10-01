from sunpy.net.attr import SimpleAttr

__all__ = ["SatelliteNumber"]


# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__


class SatelliteNumber(SimpleAttr):
    """
    The GOES Satellite Number
    """
