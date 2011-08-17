"""SDO Map subclass definitions"""
#pylint: disable=W0221

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap
from datetime import datetime
from sunpy.cm import cm
from sunpy.util import util as util

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
            'date': util.anytim(header['date-obs'][0:22]),
            'det': "AIA",
            'inst': "AIA",
            'meas': header['wavelnth'],
            'obs': "SDO",
            'name': "AIA %s" % header['wavelnth'],
            'r_sun': header['rsun_obs'],
            'cmap': cm.get_cmap(name = 'sdoaia' + str(header['wavelnth']))
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume') and header.get('instrume')[0:3] == 'AIA'
        
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
            "norm": None,
            "date": util.anytim(header['date-obs'][0:22]),
            "det": "HMI",
            "inst": "HMI",
            "meas": meas,
            "obs": "SDO",
            "name": "HMI %s" % meas,
            "r_sun": header['rsun_obs']
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an HMI image"""
        return header.get('instrume') and header.get('instrume')[0:3] == 'HMI'

