"""Hinode XRT Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

import numpy as np
from matplotlib import colors

from sunpy.map import GenericMap
from sunpy.cm import cm


__all__ = ['XRTMap']

def _lower_list(L):
    return [item.lower() for item in L]

class XRTMap(GenericMap):
    """XRT Image Map definition
    
    References
    ----------
    For a description of XRT headers
    """
    #TODO: get a link for the XRT FITS headers
    # Add in some information about the the possible filter wheel measurements
    filter_wheel1_measurements = ["Al_med", "Al_poly", "Be_med",
                                  "Be_thin", "C_poly", "Open"]
    filter_wheel2_measurements = ["Open", "Al_mesh", "Al_thick",
                                  "Be_thick", "Gband", "Ti_poly"]
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        fw1 = header.get('EC_FW1_')
        if fw1.lower() not in _lower_list(self.filter_wheel1_measurements):
            raise ValueError('Unpexpected filter wheel 1 in header.')
        fw1 = fw1.replace("_", " ")    
            
        fw2 = header.get('EC_FW2_')
        if fw2.lower() not in _lower_list(self.filter_wheel2_measurements):
            raise ValueError('Unpexpected filter wheel 2 in header.')
        fw2 = fw2.replace("_", " ")
        
        self.meta['detector'] = "XRT"
#        self.meta['instrume'] = "XRT"
        self.meta['telescop'] = "Hinode"
        
        self._name = "{0} {1}-{2}".format(self.detector, fw1, fw2)
        self._nickname = self.detector
        
        self.cmap = cm.get_cmap(name='hinodexrt')

    def _get_norm(self):
        """Returns a Normalize object to be used with XRT data"""
        # byte-scaled images have most likely already been scaled
        if self.dtype == np.uint8:
            return None

        mean = self.mean()
        std = self.std()
        
        vmin = max(0, mean - 3 * std)
        vmax = min(self.max(), mean + 3 * std)
        
        return colors.Normalize(vmin, vmax)


    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an XRT image"""
        return header.get('instrume') == 'XRT'
        
