"""Yohkoh SXT Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

from sunpy.map import Map
from sunpy.cm import cm
from sunpy.sun import constants
from matplotlib import colors
import numpy as np

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
        
        # 2012/11/07 - the SXT headers do not have a value of the distance from
        # the spacecraft to the center of the Sun.  The FITS keyword 'DSUN_OBS'
        # appears to refer to the observed diameter of the Sun.  Until such 
        # time as that is calculated and properly included in the file, we will 
        # use the value of 1 AU as a standard.
        properties['dsun']= constants.au
        
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
