"""SOHO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.sun import constants
from sunpy.cm import cm
from sunpy.time import parse_time
from matplotlib import colors
import numpy as np

class EITMap(BaseMap):
    """EIT Image Map definition"""
    @classmethod
    def get_properties(cls, header):
        """Parses EIT image header"""
        properties = BaseMap.get_properties(header)
        
        # Solar radius in arc-seconds at 1 au
        # @TODO: use sunpy.sun instead
        radius_1au = 959.644
        
        scale = header.get("cdelt1")
        
        properties.update({
            "date": parse_time(header.get('date_obs')),
            "detector": "EIT",
            "dsun": ((radius_1au / 
                      (properties['rsun_arcseconds'] * scale)) * constants.au),
            "name": "EIT %s" % header.get('wavelnth'),
            "nickname": "EIT",
            "cmap": cm.get_cmap('sohoeit%d' % header.get('wavelnth'))
        })
        return properties

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an EIT image"""
        return header.get('instrume') == 'EIT'

    def norm(self):
        """Returns a Normalize object to be used with EIT data"""
        # byte-scaled images have most likely already been scaled
        if self.dtype == np.uint8:
            return None
        
        mean = self.mean()
        std = self.std()
        
        vmin = 1
        vmax = min(self.max(), mean + 5 * std)

        return colors.LogNorm(vmin, vmax)

class LASCOMap(BaseMap):
    """LASCO Image Map definition"""
    @classmethod
    def get_properties(cls, header):
        """Parses LASCO image header"""
        properties = BaseMap.get_properties(header)
        
        datestr = "%sT%s" % (header.get('date_obs'), header.get('time_obs'))
        
        properties.update({
            "date": parse_time(datestr),
            "measurement": "white-light",
            "name": "LASCO %s" % header.get('detector'),
            "nickname": "LASCO-%s" % header.get('detector'),
            "cmap": cm.get_cmap('soholasco%s' % properties['detector'][1])
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an LASCO image"""
        return header.get('instrume') == 'LASCO'
        
class MDIMap(BaseMap):
    """MDI Image Map definition"""
    @classmethod
    def get_properties(cls, header):
        """Parses MDI image header"""
        properties = BaseMap.get_properties(header)
        
        # MDI sometimes has an "60" in seconds field
        datestr = header['date_obs']

        if datestr[17:19] == "60":
            datestr = datestr[:17] + "30" + datestr[19:]
            
        rsun = header.get('radius')
        
        # Solar radius in arc-seconds at 1 au
        # @TODO: use sunpy.sun instead
        radius_1au = 959.644
        
        # MDI images may have radius = 0.0
        if not rsun:
            dsun = constants.au
        else:
            scale = header.get("cdelt1")
            dsun = (radius_1au / (rsun * scale)) * constants.au
            
        # Determine measurement
        dpcobsr = header.get('dpc_obsr')
        meas = "magnetogram" if dpcobsr.find('Mag') != -1 else "continuum"
        
        properties.update({
            "date": parse_time(datestr),
            "detector": "MDI",
            "measurement": meas,
            "dsun": dsun,
            "name": "MDI %s" % meas,
            "nickname": "MDI"
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an MDI image"""
        return header.get('instrume') == 'MDI'
