"""SDO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
from matplotlib import colors

from sunpy.map import Map
from sunpy.cm import cm

__all__ = ['AIAMap', 'HMIMap']

class AIAMap(Map):
    """AIA Image Map definition
    
    References
    ----------
    For a description of AIA headers
    http://jsoc.stanford.edu/doc/keywords/AIA/AIA02840_A_AIA-SDO_FITS_Keyword_Documents.pdf
    """
    @classmethod
    def get_properties(cls, header):
        """Parses AIA image header"""
        properties = Map.get_properties(header)
        
        properties.update({
            "detector": "AIA",
            "instrument": "AIA",
            "observatory": "SDO",
            "nickname": "AIA",
            "cmap": cm.get_cmap('sdoaia%d' % header.get('wavelnth')),
            "processing_level": header.get('LVL_NUM')            
        })
        return properties

    def norm(self):
        """Returns a Normalize object to be used with AIA data"""
        # byte-scaled images have most likely already been scaled
        if self.dtype == np.uint8:
            return None

        mean = self.mean()
        std = self.std()
        
        vmin = max(0, mean - 3 * std)
        vmax = min(self.max(), mean + 3 * std)
        
        return colors.Normalize(vmin, vmax)
    
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume', '').startswith('AIA')
        
class HMIMap(Map):
    """HMI Image Map definition"""
    @classmethod
    def get_properties(cls, header):
        """Parses HMI image header"""
        properties = Map.get_properties(header)
        
        measurement = header['content'].split(" ")[0].lower()
        
        properties.update({
            "detector": "HMI",
            "instrument": "HMI",
            "measurement": measurement,
            "observatory": "SDO",
            "name": "HMI %s" % measurement,
            "nickname": "HMI"
        })
        return properties

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an HMI image"""
        return header.get('instrume', '').startswith('HMI') 
