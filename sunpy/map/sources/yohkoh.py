"""Yohkoh SXT Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

from sunpy.map import Map
from sunpy.cm import cm

class SXTMap(Map):
    """SXT Image Map definition
    
    References
    ----------
    For a description of SXT headers
    http://proba2.oma.be/index.html/swap/swap-analysis-manual/article/data-products?menu=23
    """
    @classmethod
    def get_properties(cls, header):
        """Parses SXT image header"""
        properties = Map.get_properties(header)
        
        wavelnth = header.get('wavelnth')
        if wavelnth == 'Al.1':
            wavelnth = 'Al01'
        
        properties.update({
            "detector": "SXT",
            "instrument": "SXT",
            "observatory": "Yohkoh",
            "name": "SXT %s" % wavelnth,
            "nickname": "SXT",
            "cmap": cm.get_cmap(name='yohkohsxt' + wavelnth.lower())
        })
        return properties 

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an SXT image"""
        return header.get('instrume') == 'SXT'
