"""STEREO Map subclass definitions

Author: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.data.map import BaseMap
from datetime import datetime
import matplotlib.colors as colors
import matplotlib.cm as cm

date_format = "%Y-%m-%dT%H:%M:%S.%f"

class EUVIMap(BaseMap):
    """EUVI Image Map definition"""
    def __new__(self, data, header):
        self.cmap = cm.gray
        self.norm = colors.Normalize(5, 1024, True)
        self.date = datetime.strptime(header['date_obs'], date_format)
        self.det = "EUVI"
        self.inst = "SECCHI"
        self.meas = header['wavelnth']
        self.obs = header['obsrvtry']
        self.name = "EUVI %s" % header['wavelnth']
        self.r_sun = header['rsun']
        
        return BaseMap.__new__(self, data, header)
        
    @classmethod
    def is_datasource_for(self, header):
        """Determines if header corresponds to an EUVI image"""
        return header['detector'] == 'EUVI'
        
class CORMap(BaseMap):
    """COR Image Map definition"""
    def __new__(self, data, header):
        self.cmap = cm.gray
        self.norm = colors.Normalize(5, 1024, True)
        self.date = datetime.strptime(header['date_obs'], date_format)
        self.det = header['detector']
        self.inst = "SECCHI"
        self.meas = header['wavelnth']
        self.obs = header['obsrvtry']
        self.name = "SECCHI %s" % header['detector']
        self.r_sun = header['rsun']
        
        return BaseMap.__new__(self, data, header)
        
    @classmethod
    def is_datasource_for(self, header):
        """Determines if header corresponds to an COR image"""
        return (header['detector'] == 'COR1') or (header['detector'] == 'COR2')
