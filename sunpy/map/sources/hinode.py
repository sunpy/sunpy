"""Hinode XRT Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

from sunpy.map import Map
from sunpy.cm import cm

class XRTMap(Map):
    """XRT Image Map definition
    
    References
    ----------
    For a description of XRT headers
    http://proba2.oma.be/index.html/swap/swap-analysis-manual/article/data-products?menu=23
    """
    @classmethod
    def get_properties(cls, header):
        """Parses XRT image header"""
        properties = Map.get_properties(header)
        
        wavelnth = header.get('wavelnth')
        if wavelnth == 'Al.1':
            wavelnth = 'Al01'
        
        properties.update({
            "detector": "XRT",
            "instrument": "XRT",
            "observatory": "Hinode",
            "name": "Hinode %s" % wavelnth,
            "nickname": "XRT",
            "cmap": cm.get_cmap(name='hinodexrt' + wavelnth.lower())
        })
        return properties 

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an XRT image"""
        return header.get('instrume') == 'XRT'
