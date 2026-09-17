from sunpy.net.attr import SimpleAttr

__all__ = ['ADAPTFileType', 'ADAPTLonType', 'ADAPTInputSource',
           'ADAPTDataAssimilation', 'ADAPTResolution', 'ADAPTVersionYear', 'ADAPTVersionMonth',
           'ADAPTEvolutionMode', 'ADAPTHelioData', 'ADAPTRealizations', 'ADAPTMagData']

# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__


class ADAPTFileType(SimpleAttr):
    """
    ADAPT file type: Public.
    """

class ADAPTLonType(SimpleAttr):
    """
    ADAPT longitude type: Carrington Fixed, Central Meridian, East Limb.
    """

class ADAPTInputSource(SimpleAttr):
    """
    ADAPT input source: All, KPVT, VSM, GONG, HMI, MDI, MWO.
    """

class ADAPTDataAssimilation(SimpleAttr):
    """
    ADAPT data assimilation: WH, enLS, enkf, enLAKF.
    """

class ADAPTResolution(SimpleAttr):
    """
    ADAPT model spatial resolution: 1.0 deg, 0.2 deg.
    """

class ADAPTVersionYear(SimpleAttr):
    """
    ADAPT code version year.
    """

class ADAPTVersionMonth(SimpleAttr):
    """
    ADAPT code version month.
    """

class ADAPTEvolutionMode(SimpleAttr):
    """
    ADAPT evolution mode: Data assimilation step, Intermediate step, Forecast step.
    """

class ADAPTHelioData(SimpleAttr):
    """
    ADAPT helioseismic data: Not added or no data, Far-side, Emergence, Both emergence & far-side.
    """

class ADAPTRealizations(SimpleAttr):
    """
    ADAPT realizations: number of model/file realizations, e.g., 16 -> "016"
    """


class ADAPTMagData(SimpleAttr):
    """
    ADAPT magnetic data: Not added or no data, Mag-los, Mag-vector, Mag- both los & vector, Mag- polar avg obs, Mag- los & polar, Mag- vector & polar, Mag- both los and vector & polar.
    """
