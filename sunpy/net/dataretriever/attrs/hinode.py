from sunpy.net.attr import SimpleAttr

__all__ = ['SOTDetector']


# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__


class SOTDetector(SimpleAttr):
    """
    Solar Optical Telescope (SOT) instrument, e.g., 'FG' or 'SP'.
    """
    pass
