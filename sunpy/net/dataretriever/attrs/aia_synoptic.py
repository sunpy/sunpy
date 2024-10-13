from sunpy.net.attr import SimpleAttr

__all__ = ['ADAPTFileType', 'ADAPTLonType', 'ADAPTInputSource',
           'ADAPTDataAssimilation', 'ADAPTResolution', 'ADAPTVersionYear', 'ADAPTVersionMonth',
           'ADAPTEvolutionMode', 'ADAPTHelioData', 'ADAPTRealizations', 'ADAPTMagData']

# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__


class AIASynopticData(SimpleAttr):
    """
    Custom attribute to indicate the use of low-resolution synoptic AIA data.
    Defaults to True if no value is provided.
    """

    def __init__(self, value=True):
        super().__init__(value)

    def __repr__(self):
        return f"AIASynopticData({self.value})"
