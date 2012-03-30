"""PROBA2 Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101

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

    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        properties = BaseMap.get_properties(header)
        properties.update({
            'date': parse_time(header.get('date-obs')),
            'det': "SWAP",
            'inst': "SWAP",
            'meas': header.get('wavelnth'),
            'obs': "PROBA2",
            'dsun': header.get('dsun_obs'),
            'name': "SWAP %s" % header.get('wavelnth'),
            'cmap': cm.get_cmap(name='sdoaia171'),
            'exptime': header.get('exptime')
        })
        return properties
    
    def norm(self):
        """Returns a Normalize object to be used with SWAP data"""
        mean = self.mean()
        std = self.std()
        
        vmin = max(0, mean - 3 * std)
        vmax = min(self.max(), mean + 3 * std)

        return colors.Normalize(vmin, vmax)
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an SWAP image"""
        return header.get('instrume') and header.get('instrume') == 'SWAP'
