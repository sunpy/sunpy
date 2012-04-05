"""PROBA2 Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121,W0613

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.cm import cm
from sunpy.time import parse_time
from matplotlib import colors

class SWAPMap(BaseMap):
    """SWAP Image Map definition
    
    Reference
    ---------
    For a description of SWAP headers
    http://proba2.oma.be/index.html/swap/swap-analysis-manual/article/data-products?menu=23
    """
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
    
    def __init__(self, data, header):
        BaseMap.__init__(self, header)
        self.det = "SWAP"
        self.inst = "SWAP"
        self.obs = "PROBA2"
        self.name = "SWAP %s" % header.get('wavelnth')
        self.cmap = cm.get_cmap(name='sdoaia171')

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an SWAP image"""
        return header.get('instrume') == 'SWAP'
