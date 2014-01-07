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

def _dsunAtSoho(date, rad_d, rad_1au = None):
    """Determines the distance to the Sun from SOhO following
    d_{\sun,Object} =
            D_{\sun\earth} \frac{\tan(radius_{1au}[rad])}{\tan(radius_{d}[rad])}
    though tan x ~ x for x << 1
    d_{\sun,Object} =
            D_{\sun\eart} \frac{radius_{1au}[rad]}{radius_{d}[rad]}
    since radius_{1au} and radius_{d} are dividing each other we can use [arcsec]
    instead. 

    ---
    TODO: Does this apply just to observations on the same Earth-Sun line?
    If not it can be moved outside here.
    """
    if not rad_1au:
        rad_1au = sun.solar_semidiameter_angular_size(date)
    dsun = sun.sunearth_distance(date) * constants.au * (rad_1au / rad_d)
    #return scalar value not astropy.quantity
    return dsun.value


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
        dsun = _dsunAtSoho(self.date, self.rsun_arcseconds)
        self.meta['dsun_obs'] = dsun

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EIT image"""
        return header.get('instrume') == 'EIT'

    def _get_norm(self):
        """Returns a Normalize object to be used with EIT data"""
        # byte-scaled images have most likely already been scaled
        # THIS WARNING IS KNOWN TO APPLY TO 0.3 code only.
        # NOT TESTED in sunpy 0.4 when the glymur library
        # is used instead of pyopenjpeg.  It seems that EIT JP2 files read by
        # pyopenjpeg and openjpeg using the j2k_to_image command, returns 
        # np.float32 arrays.  For comparison, AIA JP2 files read the same way
        # return np.uint8 arrays.  EIT JP2 files have already been
        # byte-scaled when they are created by the Helioviewer Project.
        # SunPy 0.3 and lower code assumes that if datatype of the data array
        # was np.uint8 then the image was highly likely to be byte-scaled
        # Since the data datatype was in fact a np.float32, then the byte-scaling
        # was never picked up.
        if self.data.dtype == np.float32:
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
        if 'date_obs' not in self.meta:
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
            scale = self.meta.get('cdelt1')
            self.meta['dsun_obs'] = _dsunAtSoho(self.date, radius * scale)
        
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an MDI image"""
        return header.get('instrume') == 'MDI'

