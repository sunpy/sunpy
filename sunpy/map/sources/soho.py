"""SOHO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
from matplotlib import colors

from sunpy.map import GenericMap
from sunpy.sun import constants
from sunpy.sun import sun
from sunpy.cm import cm
from sunpy.time import parse_time, is_time

__all__ = ['EITMap', 'LASCOMap', 'MDIMap']

class EITMap(GenericMap):
    """EIT Image Map definition"""
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        # Fill in some missing info
        self.meta['detector'] = "EIT"
        self._fix_dsun()
        
        self._name = self.detector + " " + str(self.measurement)
        self._nickname = self.detector
        
        self.cmap = cm.get_cmap('sohoeit%d' % self.wavelength)
    
    @property
    def rsun_arcseconds(self):
        return self.meta['solar_r'] * self.meta['cdelt1']
        
    def _fix_dsun(self):
        dsun = constants.au * sun.angular_size(self.date) / (self.rsun_arcseconds)
        self.meta['dsun_obs'] = dsun

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
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

class LASCOMap(GenericMap):
    """LASCO Image Map definition"""
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        # Fill in some missing or broken info
        self._fix_date()
        
        self._name = self.instrument + " " + self.detector
        self._nickname = self.instrument + "-" + self.detector
                
        self.cmap = cm.get_cmap('soholasco%s' % self.detector[1])
        
    @property
    def measurement(self):
        # TODO: This needs to do more than white-light.  Should give B, pB, etc.
        return "white-light"
    
    def _fix_date(self):
        datestr = "%sT%s" % (self.meta.get('date-obs',self.meta.get('date_obs')),
                     self.meta.get('time-obs',self.meta.get('time_obs')))

        # If non-standard Keyword is present, correct it too, for compatibility.
        if 'date_obs' in self.meta:
            self.meta['date_obs'] = self.meta['date-obs']
            
        
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an LASCO image"""
        return header.get('instrume') == 'LASCO'
        
class MDIMap(GenericMap):
    """MDI Image Map definition"""
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        # Fill in some missing or broken info
        self.meta['detector'] = "MDI"
        self._fix_dsun()
        
        self._name = self.observatory + " " + self.detector
        self._nickname = self.detector + " " + self.measurement
        
    @property
    def measurement(self):
        # TODO: This needs to do more than white-light.  Should give B, pB, etc.
        return "magnetogram" if self.meta['dpc_obsr'].find('Mag') != -1 else "continuum"
        
    def _fix_dsun(self):
        # Solar radius in arc-seconds at 1 au
        # previous value radius_1au = 959.644
        # radius = constants.average_angular_size
        radius = self.meta.get('radius')
        
        if not radius:
#            radius = sun.angular_size(self.date)
            self.meta['dsun_obs'] = constants.au
        else:
            self.meta['dsun_obs'] = constants.au * radius / (self.rsun_arcseconds)
        
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an MDI image"""
        return header.get('instrume') == 'MDI'

