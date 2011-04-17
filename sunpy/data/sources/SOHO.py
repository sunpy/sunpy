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
        return BaseMap.__new__(self, data, header)
        
    @classmethod
    def get_properties(self, header):
        """Returns the default and normalized values to use for the Map"""
        date_format = "%Y-%m-%dT%H:%M:%S.%fZ"
        
        return {
            "cmap": cm.gray,
            "norm": colors.Normalize(5, 1024, True),
            "date": datetime.strptime(header['date_obs'], date_format),
            "det": "EIT",
            "inst": "EIT",
            "meas": header['wavelnth'],
            "obs": "SOHO",
            "name": "EIT %s" % header['wavelnth'],
            "r_sun": header['solar_r']
        }
        
    @classmethod
    def is_datasource_for(self, header):
        """Determines if header corresponds to an EIT image"""
        return header['instrume'] == 'EIT'

class LASCOMap(BaseMap):
    """LASCO Image Map definition"""
    def __new__(self, data, header):
        return BaseMap.__new__(self, data, header)
        
    @classmethod
    def get_properties(self, header):
        """Returns the default and normalized values to use for the Map"""
        datestr = "%sT%s" % (header['date_obs'], header['time_obs'])

        return {
            "cmap": cm.gray,
            "norm": colors.Normalize(5, 1024, True),
            "date": datetime.strptime(datestr, "%Y/%m/%dT%H:%M:%S.%f"),
            "det": header['detector'],
            "inst": "LASCO",
            "meas": header['wavelnth'],
            "obs": "SOHO",
            "name": "LASCO %s" % header['detector'],
            "r_sun": None
        }
        
    @classmethod
    def is_datasource_for(self, header):
        """Determines if header corresponds to an LASCO image"""
        return header['instrume'] == 'LASCO'
        
class MDIMap(BaseMap):
    """MDI Image Map definition"""
    def __new__(self, data, header):
        return BaseMap.__new__(self, data, header)
        
    @classmethod
    def get_properties(self, header):
        """Returns the default and normalized values to use for the Map"""
        # MDI sometimes has an "60" in seconds field
        datestr = header['date_obs']

        if datestr[17:19] == "60":
            datestr = datestr[:17] + "30" + datestr[19:]
        
        # Determine measurement
        dpcobsr = header['dpc_obsr']
        meas = "magnetogram" if dpcobsr.find('Mag') != -1 else "continuum"
        
        return {
            "cmap": cm.gray,
            "norm": colors.Normalize(5, 1024, True),
            "date": datetime.strptime(datestr, "%Y-%m-%dT%H:%M:%S.%fZ"),
            "det": "MDI",
            "inst": "MDI",
            "meas": meas,
            "obs": "SOHO",
            "name": "MDI %s" % meas,
            "r_sun": header['r_sun']
        }
        
    @classmethod
    def is_datasource_for(self, header):
        """Determines if header corresponds to an MDI image"""
        return header['instrume'] == 'MDI'

