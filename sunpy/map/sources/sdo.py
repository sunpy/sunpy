"""SDO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy import visualization

from sunpy.map import GenericMap
from sunpy.cm import cm

__all__ = ['AIAMap', 'HMIMap']

class AIAMap(GenericMap):
    """AIA Image Map.

    The Atmospheric Imaging Assembly is a set of four telescopes that employ
    normal-incidence, multi-layer coated optics to provide narrow-band imaging
    of the Sun. It provides high resolution full-disk images of the corona and
    transition region up to 0.5 solar radii above the solar limb with 1.5
    arcsecond angular resolution and 12-second temporal resolution. It observes
    the Sun in the following seven extreme ultraviolet bandpasses: 94 A
    (Fe XVIII), 131 A (Fe VIII, XXI), 171 A (Fe IX), 193 A (Fe XII, XXIV),
    211 A (Fe XIV), 304 A (He II), 335 A (Fe XVI). One telescope observes
    in the visible 1600 A (C IV) and the nearby continuun (1700 A).

    References
    ----------
    * `SDO Mission Page <http://sdo.gsfc.nasa.gov>`_
    * `Instrument Page <http://aia.lmsal.com>`_
    * `Fits Header keywords <http://jsoc.stanford.edu/doc/keywords/AIA/AIA02840_A_AIA-SDO_FITS_Keyword_Documents.pdf>`_
    * `Analysis Guide <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/>`_
    * `Instrument Paper <http://link.springer.com/article/10.1007%2Fs11207-011-9776-8>`_
    * `wavelengths and temperature response reference <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/figures/aia_tel_resp.png>`_
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        # Fill in some missing info
        self.meta['detector'] = "AIA"
        self._nickname = self.detector
        self.plot_settings['cmap'] = cm.get_cmap(self._get_cmap_name())
        self.plot_settings['norm'] = ImageNormalize(stretch=visualization.AsinhStretch(0.01))

    @property
    def observatory(self):
        """
        Returns the observatory.
        """
        return self.meta['telescop'].split('/')[0]

    @property
    def processing_level(self):
        """
        Returns the FITS processing level.
        """
        return self.meta['lvl_num']

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume', '').startswith('AIA')

class HMIMap(GenericMap):
    """HMI Image Map.

    HMI consists of a refracting telescope, a polarization selector,
    an image stabilization system, a narrow band tunable filter
    and two 4096 pixel CCD cameras. It observes the full solar disk in the Fe I
    absorption line at 6173 Angstrom with a resolution of 1 arc-second.
    HMI takes images in a sequence of tuning and polarizations at a 4-second
    cadence for each camera. One camera is dedicated to a 45 s Doppler and
    line-of-sight field sequence while the other to a 90 s vector field
    sequence.

    References
    ----------
    * `SDO Mission Page <http://sdo.gsfc.nasa.gov>`_
    * `Instrument Page <http://hmi.stanford.edu>`_
    * `Analysis Guide <http://hmi.stanford.edu/doc/magnetic/guide.pdf>`_
    """
    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        self.meta['detector'] = "HMI"
#        self.meta['instrme'] = "HMI"
#        self.meta['obsrvtry'] = "SDO"
        self.meta['waveunit'] = 'Angstrom'
        self._nickname = self.detector

    @property
    def measurement(self):
        """
        Returns the measurement type.
        """
        return self.meta['content'].split(" ")[0].lower()

    @property
    def observatory(self):
        """
        Returns the observatory.
        """
        return self.meta['telescop'].split('/')[0]

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an HMI image"""
        return header.get('instrume', '').startswith('HMI')
