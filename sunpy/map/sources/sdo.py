"""SDO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.cm import cm
from sunpy.util import util as util
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

    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        # Note: Trailing "Z" in date was dropped on 2010/12/07        
        properties = BaseMap.get_properties()
        properties.update({
            'date': util.anytim(header.get('date-obs')),
            'det': "AIA",
            'inst': "AIA",
            'meas': header.get('wavelnth'),
            'obs': "SDO",
            'name': "AIA %s" % header.get('wavelnth'),
            'cmap': cm.get_cmap(name = 'sdoaia' + str(header.get('wavelnth'))),
            'exptime': header.get('exptime')
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume') and header.get('instrume')[0:3] == 'AIA'
    
    def norm(self):
        """Returns a Normalize object to be used with AIA data"""
        mean = self.mean()
        std = self.std()
        
        vmin = max(0, mean - 3 * std)
        vmax = min(self.max(), mean + 3 * std)
        
        # 8-bit images are probably from Helioviewer and are already scaled
        vmax = max(255, vmax)
        
        return colors.Normalize(vmin, vmax)
        
class HMIMap(BaseMap):
    """HMI Image Map definition"""
    def __new__(cls, data, header):        
        return BaseMap.__new__(cls, data)
        
    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        # Note: Trailing "Z" in date was dropped on 2010/12/07    
        meas = header['content'].split(" ")[0].lower()
        
        properties = BaseMap.get_properties()
        properties.update({
            "date": util.anytim(header.get('date-obs')),
            "det": "HMI",
            "inst": "HMI",
            "meas": meas,
            "obs": "SDO",
            "name": "HMI %s" % meas,
            "r_sun": header.get('rsun_obs'),
            "exptime": header.get('exptime')
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an HMI image"""
        return header.get('instrume') and header.get('instrume')[0:3] == 'HMI'

