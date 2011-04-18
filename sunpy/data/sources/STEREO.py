"""STEREO Map subclass definitions

Author: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.data.BaseMap import BaseMap
from datetime import datetime

class EUVIMap(BaseMap):
    """EUVI Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data, header)
        
    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        date_format = "%Y-%m-%dT%H:%M:%S.%f"

        properties = BaseMap.get_properties()
        properties.update({
            "date": datetime.strptime(header['date_obs'], date_format),
            "det": "EUVI",
            "inst": "SECCHI",
            "meas": header['wavelnth'],
            "obs": header['obsrvtry'],
            "name": "EUVI %s" % header['wavelnth'],
            "r_sun": header['rsun']
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an EUVI image"""
        return header['detector'] == 'EUVI'
        
class CORMap(BaseMap):
    """COR Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data, header)
        
    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        properties = BaseMap.get_properties()
        date_format = "%Y-%m-%dT%H:%M:%S.%f"

        properties.update({
            "date": datetime.strptime(header['date_obs'], date_format),
            "det": header['detector'],
            "inst": "SECCHI",
            "meas": header['wavelnth'],
            "obs": header['obsrvtry'],
            "name": "SECCHI %s" % header['detector'],
            "r_sun": header['rsun']
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an COR image"""
        return (header['detector'] == 'COR1') or (header['detector'] == 'COR2')
