"""SDO Map subclass definitions

Author: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.data.map import BaseMap
from datetime import datetime
import matplotlib.colors as colors
import matplotlib.cm as cm

# Note: Trailing "Z" in date was dropped on 2010/12/07
date_format = "%Y-%m-%dT%H:%M:%S.%f"

class AIAMap(BaseMap):
    """AIA Image Map definition"""
    def __new__(self, data, header):
        self.cmap = cm.gray
        self.norm = colors.Normalize(5, 1024, True)
        self.date = datetime.strptime(header['date-obs'][0:22], date_format)
        self.det = "AIA"
        self.inst = "AIA"
        self.meas = header['wavelnth']
        self.obs = "SDO"
        self.name = "AIA %s" % header['wavelnth']
        self.r_sun = header['rsun_obs']
        
        return BaseMap.__new__(self, data, header)
        
    @classmethod
    def is_datasource_for(self, header):
        """Determines if header corresponds to an AIA image"""
        return header['telescop'] == 'SDO/AIA'
        
class AIAMap(BaseMap):
    """AIA Image Map definition"""
    def __new__(self, data, header):
        meas = header['content'].split(" ")[0].lower()
        self.cmap = cm.gray
        self.norm = None
        self.date = datetime.strptime(header['date-obs'][0:22], date_format)
        self.det = "HMI"
        self.inst = "HMI"
        self.meas = meas
        self.obs = "SDO"
        self.name = "HMI %s" % meas
        self.r_sun = header['rsun_obs']
        
        return BaseMap.__new__(self, data, header)
        
    @classmethod
    def is_datasource_for(self, header):
        """Determines if header corresponds to an HMI image"""
        return header['instrumu'][0:3] == 'HMI'

