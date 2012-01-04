"""SOHO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.cm import cm
from sunpy.util import util
from matplotlib import colors

class EITMap(BaseMap):
    """EIT Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
        
    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        #date_format = "%Y-%m-%dT%H:%M:%S.%fZ"
        
        properties = BaseMap.get_properties()
        properties.update({
            "date": util.anytim(header.get('date_obs')),
            "det": "EIT",
            "inst": "EIT",
            "meas": header.get('wavelnth'),
            "obs": "SOHO",
            "name": "EIT %s" % header.get('wavelnth'),
            "exptime": header.get('exptime'),
            'cmap': cm.get_cmap(name='sohoeit' + str(header.get('wavelnth'))),
            "r_sun": header.get('solar_r')
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

        properties = BaseMap.get_properties()
        properties.update({
            "date": util.anytim(datestr),
            "det": header.get('detector'),
            "inst": "LASCO",
            "meas": header.get('wavelnth'),
            "obs": "SOHO",
            "name": "LASCO %s" % header.get('detector'),
            "exptime": header.get("exptime"),
            "r_sun": None
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
        
        # Determine measurement
        dpcobsr = header.get('dpc_obsr')
        meas = "magnetogram" if dpcobsr.find('Mag') != -1 else "continuum"
        
        properties = BaseMap.get_properties()        
        properties.update({
            "date": util.anytim(datestr),
            "det": "MDI",
            "inst": "MDI",
            "meas": meas,
            "obs": "SOHO",
            "name": "MDI %s" % meas,
            "exptime": header.get('exptime'),
            "r_sun": header.get('radius')
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an MDI image"""
        return header.get('instrume') == 'MDI'
