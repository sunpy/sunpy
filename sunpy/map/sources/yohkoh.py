"""Yohkoh SXT Map subclass definitions"""

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

import numpy as np

from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch
from sunpy.sun import constants

__all__ = ['SXTMap']


class SXTMap(GenericMap):
    """Yohkoh SXT Image Map

    The Yohkoh Soft X-ray Telescope (SXT) the full solar disk
    (42 x 42 arcminutes)in the 0.25 - 4.0 keV range.
    It consists of a glancing incidence mirror and a CCD sensor and
    used thin metallic filters to acquire images in restricted
    portions of its energy range. SXT could resolve features down to 2.5
    arcseconds. Information about the temperature and density of the plasma
    emitting the observed x-rays was obtained by comparing images acquired with
    the different filters. Images could be obtained every 2 to 8 seconds.
    Smaller images with a single filter could be obtained as frequently as
    once every 0.5 seconds.

    Yohkoh was launched on 30 August 1991 and ceased operations on
    14 December 2001.

    References
    ----------
    * `Yohkoh Mission Page <http://solar.physics.montana.edu/sxt/>`_
    * `Fits header reference <http://proba2.oma.be/data/SWAP/level0>`_
    * `Yohkoh Analysis Guide <http://ylstone.physics.montana.edu/ylegacy/yag.html>`_
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        self.meta['detector'] = "SXT"
        self.meta['telescop'] = "Yohkoh"
        self.plot_settings['cmap'] = 'yohkohsxt' + self.measurement[0:2].lower()
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, PowerStretch(0.5)), clip=False)

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
    def measurement(self):
        """
        Returns the type of data observed.
        """
        s = self.meta.get('wavelnth', '')
        if s == 'Al.1':
            s = 'Al01'
        elif s.lower() == 'open':
            s = 'white-light'
        return s

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an SXT image"""
        return header.get('instrume') == 'SXT'
