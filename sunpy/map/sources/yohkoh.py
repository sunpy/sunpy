"""Yohkoh SXT Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

import numpy as np
from matplotlib import colors

from sunpy.map import GenericMap
from sunpy.cm import cm
from sunpy.sun import constants

__all__ = ['SXTMap']

class SXTMap(GenericMap):
    """SXT Image Map definition
    
    References
    ----------
    For a description of SXT headers
    http://proba2.oma.be/index.html/swap/swap-analysis-manual/article/data-products?menu=23
    """
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        self.meta['detector'] = "SXT"
        self.meta['telescop'] = "Yohkoh"
        
        self._name = self.observatory + " " + self.wavelength_string
    
        self.cmap = cm.get_cmap(name='yohkohsxt' + self.wavelength_string[0:2].lower())
    
        # 2012/12/19 - the SXT headers do not have a value of the distance from
        # the spacecraft to the center of the Sun.  The FITS keyword 'DSUN_OBS'
        # appears to refer to the observed diameter of the Sun.  Until such 
        # time as that is calculated and properly included in the file, we will 
        # use simple trigonometry to calculate the distance of the center of 
        # the Sun from the spacecraft.  Note that the small angle approximation
        # is used, and the solar radius stored in SXT FITS files is in arcseconds.
        self.meta['dsun_apparent'] = constants.au
        if 'solar_r' in self.meta:
            self.meta['dsun_apparent'] = constants.radius/(np.deg2rad(self.meta['solar_r']/3600.0))
   
    @property
    def dsun(self):
        """ For Yohkoh Maps, dsun_obs is not always defined. Uses approximation
        defined above it is not defined."""
        return self.meta.get('dsun_obs', self.meta['dsun_apparent'])
    
    @property
    def wavelength_string(self):
        s = self.meta.get('wavelnth', '')
        if s == 'Al.1':
            s = 'Al01' 
        elif s.lower() ==  'open':
            s = 'white light'
        return s

    def _get_norm(self):
        """Returns a Normalize object to be used with SXT data"""
        # byte-scaled images have most likely already been scaled
        if self.dtype == np.uint8:
            return None

        mean = self.mean()
        std = self.std()
        
        vmin = max(0, mean - 3 * std)
        vmax = min(self.max(), mean + 3 * std)
        
        return colors.Normalize(vmin, vmax)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an SXT image"""
        return header.get('instrume') == 'SXT'
