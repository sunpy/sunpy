"""PROBA2 Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map import Map
from sunpy.cm import cm

__all__ = ['SWAPMap']

class SWAPMap(Map):
    """SWAP Image Map definition
    
    References
    ----------
    For a description of SWAP headers
    http://proba2.oma.be/index.html/swap/swap-analysis-manual/article/data-products?menu=23
    """
    @classmethod
    def get_properties(cls, header):
        """Parses SWAP image header"""
        properties = Map.get_properties(header)
        
        properties.update({
            "detector": "SWAP",
            "instrument": "SWAP",
            "observatory": "PROBA2",
            "name": "SWAP %s" % header.get('wavelnth'),
            "nickname": "SWAP",
            "cmap": cm.get_cmap(name='sdoaia171')
        })
        return properties

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an SWAP image"""
        return header.get('instrume') == 'SWAP'
