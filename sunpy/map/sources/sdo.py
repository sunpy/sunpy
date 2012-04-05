"""SDO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121,W0613

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.cm import cm
from matplotlib import colors

class AIAMap(BaseMap):
    """AIA Image Map definition
    
    Reference
    ---------
    For a description of AIA headers
    http://jsoc.stanford.edu/doc/keywords/AIA/AIA02840_A_AIA-SDO_FITS_Keyword_Documents.pdf
    """
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)

    def __init__(self, data, header):
        BaseMap.__init__(self, header)
        
        self.det = "AIA"
        self.inst = "AIA"
        self.obs = "SDO"
        self.cmap = cm.get_cmap('sdoaia%d' % header.get('wavelnth'))

    def norm(self):
        """Returns a Normalize object to be used with AIA data"""
        mean = self.mean()
        std = self.std()
        
        vmin = max(0, mean - 3 * std)
        vmax = min(self.max(), mean + 3 * std)
        
        # 8-bit images are probably from Helioviewer and are already scaled
        vmax = max(255, vmax)
        
        return colors.Normalize(vmin, vmax)
    
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume', '').startswith('AIA')
        
class HMIMap(BaseMap):
    """HMI Image Map definition"""
    def __new__(cls, data, header):        
        return BaseMap.__new__(cls, data)
    
    def __init__(self, data, header):
        BaseMap.__init__(self, header)
        
        self.det = "HMI"
        self.inst = "HMI"
        self.meas = header['content'].split(" ")[0].lower()
        self.obs = "SDO"
        self.name = "HMI %s" % self.meas 

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an HMI image"""
        return header.get('instrume', '').startswith('HMI') 
