"""SOHO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.sun import constants
from sunpy.cm import cm
from sunpy.time import parse_time
from matplotlib import colors

class EITMap(BaseMap):
    """EIT Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
        
    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        #date_format = "%Y-%m-%dT%H:%M:%S.%fZ"
        
        # Solar radius in arc-seconds at 1 au
        # @TODO: use sunpy.sun instead
        radius_1au = 959.644
        
        rsun = header.get('solar_r')
        scale = header.get("cdelt1")
        dsun = (radius_1au / (rsun * scale)) * constants.au
        
        properties = BaseMap.get_properties(header)
        properties.update({
            "date": parse_time(header.get('date_obs')),
            "det": "EIT",
            "inst": "EIT",
            "meas": header.get('wavelnth'),
            "obs": "SOHO",
            "dsun": dsun,            
            "name": "EIT %s" % header.get('wavelnth'),
            "exptime": header.get('exptime'),
            'cmap': cm.get_cmap(name='sohoeit' + str(header.get('wavelnth'))),
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an EIT image"""
        return header.get('instrume') == 'EIT'
    
    
    def norm(self):
        """Returns a Normalize object to be used with AIA data"""
        mean = self.mean()
        std = self.std()
        
        vmin = 1
        vmax = min(self.max(), mean + 5 * std)
        
        # 8-bit images are probably from Helioviewer and are already scaled
        vmax = max(255, vmax)
        
        return colors.LogNorm(vmin, vmax)

class LASCOMap(BaseMap):
    """LASCO Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
        
    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        datestr = "%sT%s" % (header['date_obs'], header['time_obs'])
        
        # @TODO: Verify
        # No rsun/dsun information in LASCO images?
        dsun = constants.au

        properties = BaseMap.get_properties(header)
        properties.update({
            "date": parse_time(datestr),
            "det": header.get('detector'),
            "inst": "LASCO",
            "meas": "white-light",
            "obs": "SOHO",
            "dsun": dsun,
            "name": "LASCO %s" % header.get('detector'),
            "exptime": header.get("exptime"),
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an LASCO image"""
        return header.get('instrume') == 'LASCO'
        
class MDIMap(BaseMap):
    """MDI Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
        
    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
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
        
        properties = BaseMap.get_properties(header)        
        properties.update({
            "date": parse_time(datestr),
            "det": "MDI",
            "inst": "MDI",
            "meas": meas,
            "obs": "SOHO",
            "dsun": dsun,
            "name": "MDI %s" % meas,
            "exptime": header.get('exptime')
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an MDI image"""
        return header.get('instrume') == 'MDI'
