"""SDO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.cm import cm
from matplotlib import colors
import numpy as np

class AIAMap(BaseMap):
    """AIA Image Map definition
    
    References
    ----------
    For a description of AIA headers
    http://jsoc.stanford.edu/doc/keywords/AIA/AIA02840_A_AIA-SDO_FITS_Keyword_Documents.pdf
    """
    @classmethod
    def get_properties(cls, header):
        """Parses AIA image header"""
        properties = BaseMap.get_properties(header)
        
        properties.update({
            "detector": "AIA",
            "instrument": "AIA",
            "observatory": "SDO",
            "nickname": "AIA",
            "cmap": cm.get_cmap('sdoaia%d' % header.get('wavelnth'))
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
        
class HMIMap(BaseMap):
    """HMI Image Map definition"""
    @classmethod
    def get_properties(cls, header):
        """Parses HMI image header"""
        properties = BaseMap.get_properties(header)
        
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
