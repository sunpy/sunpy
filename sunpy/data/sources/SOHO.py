"""SOHO Map subclass definitions

Author: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.data.map import BaseMap
from datetime import datetime
import matplotlib.colors as colors
import matplotlib.cm as cm

date_format = "%Y-%m-%dT%H:%M:%S.%f"

class EITMap(BaseMap):
    """EIT Image Map definition"""
    def __new__(self, data, header):
        date_format = "%Y-%m-%dT%H:%M:%S.%fZ"
        
        self.cmap = cm.gray
        self.norm = colors.Normalize(5, 1024, True)
        self.date = datetime.strptime(header['date_obs'], date_format)
        self.det = "EIT"
        self.inst = "EIT"
        self.meas = header['wavelnth']
        self.obs = "SOHO"
        self.name = "EIT %s" % header['wavelnth']
        self.r_sun = header['solar_r']
        
        return BaseMap.__new__(self, data, header)
        
    @classmethod
    def is_datasource_for(self, header):
        """Determines if header corresponds to an EIT image"""
        return header['instrume'] == 'EIT'

class LASCOMap(BaseMap):
    """LASCO Image Map definition"""
    def __new__(self, data, header):
        datestr = "%sT%s" % (header['date_obs'], header['time_obs'])
        
        self.cmap = cm.gray
        self.norm = colors.Normalize(5, 1024, True)
        self.date = datetime.strptime(datestr, "%Y/%m/%dT%H:%M:%S.%f")
        self.det = header['detector']
        self.inst = "LASCO"
        self.meas = header['wavelnth']
        self.obs = "SOHO"
        self.name = "LASCO %s" % header['detector']
        self.r_sun = None
        
        return BaseMap.__new__(self, data, header)
        
    @classmethod
    def is_datasource_for(self, header):
        """Determines if header corresponds to an LASCO image"""
        return header['instrume'] == 'LASCO'
        
class MDIMap(BaseMap):
    """MDI Image Map definition"""
    def __new__(self, data, header):
        # MDI sometimes has an "60" in seconds field
        datestr = header['date_obs']

        if datestr[17:19] == "60":
            datestr = datestr[:17] + "30" + datestr[19:]
        
        # Determine measurement
        dpcobsr = header['dpc_obsr']
        meas = "magnetogram" if dpcobsr.find('Mag') != -1 else "continuum"
        
        self.cmap = cm.gray
        self.norm = colors.Normalize(5, 1024, True)
        self.date = datetime.strptime(datestr, "%Y-%m-%dT%H:%M:%S.%fZ"),
        self.det = "MDI"
        self.inst = "MDI"
        self.meas = meas
        self.obs = "SOHO"
        self.name = "MDI %s" % meas
        self.r_sun = header['r_sun']
        
        return BaseMap.__new__(self, data, header)
        
    @classmethod
    def is_datasource_for(self, header):
        """Determines if header corresponds to an MDI image"""
        return header['instrume'] == 'MDI'

