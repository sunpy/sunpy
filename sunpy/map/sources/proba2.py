"""PROBA2 Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map import GenericMap
from sunpy.cm import cm

__all__ = ['SWAPMap']

class SWAPMap(GenericMap):
    """SWAP Image Map definition
    
    References
    ----------
    For a description of SWAP headers
    http://proba2.oma.be/index.html/swap/swap-analysis-manual/article/data-products?menu=23
    """
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)

		# It needs to be verified that these must actually be set and are not
		# already in the header.
        self.meta['detector'] = "SWAP"
#        self.meta['instrme'] = "SWAP"
        self.meta['obsrvtry'] = "PROBA2"
        
        self._name = self.detector + " " + str(self.measurement)
        self._nickname = self.detector
        
        self.cmap = cm.get_cmap(name='sdoaia171')
            

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an SWAP image"""
        return header.get('instrume') == 'SWAP'
