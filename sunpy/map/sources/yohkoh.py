"""Yohkoh SXT Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

import numpy as np
from matplotlib import colors

from sunpy.map import Map
from sunpy.cm import cm
from sunpy.sun import constants

__all__ = ['SXTMap']

class SXTMap(Map):
    """SXT Image Map definition
    
    References
    ----------
    For a description of SXT headers
    http://proba2.oma.be/index.html/swap/swap-analysis-manual/article/data-products?menu=23
    """
    @classmethod
    def get_properties(cls, header):
        """Parses SXT image header"""
        properties = Map.get_properties(header)
        
        # 2012/12/19 - the SXT headers do not have a value of the distance from
        # the spacecraft to the center of the Sun.  The FITS keyword 'DSUN_OBS'
        # appears to refer to the observed diameter of the Sun.  Until such 
        # time as that is calculated and properly included in the file, we will 
        # use simple trigonometry to calculate the distance of the center of 
        # the Sun from the spacecraft.  Note that the small angle approximation
        # is used, and the solar radius stored in SXT FITS files is in arcseconds.
        properties['dsun']= constants.au
        yohkoh_solar_r = header.get('solar_r', None)
        if yohkoh_solar_r == None:
            properties['dsun']= constants.au
        else:
            properties['dsun'] = constants.radius/(np.deg2rad(yohkoh_solar_r/3600.0))
        
        wavelnth = header.get('wavelnth')
        if wavelnth == 'Al.1':
            wavelnth = 'Al01'
        if wavelnth.lower() == 'open':
            wavelnth = 'white light'

        properties.update({
            "detector": "SXT",
            "instrument": "SXT",
            "observatory": "Yohkoh",
            "name": "SXT %s" % wavelnth,
            "nickname": "SXT",
            "cmap": cm.get_cmap(name='yohkohsxt' + wavelnth[0:2].lower())
        })
        return properties 

    def norm(self):
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
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an SXT image"""
        return header.get('instrume') == 'SXT'
