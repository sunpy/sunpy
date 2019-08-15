"""SUVI Map subclass definitions"""
from __future__ import absolute_import, print_function, division
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

import matplotlib.pyplot as plt

from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import AsinhStretch

from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ['SUVIMap']


class SUVIMap(GenericMap):
    """SUVI Image Map.

    The Solar Ultraviolet Imager (SUVI) is a normal-incidence Cassegrain telescope
    that shares considerable design heritage with the Atmospheric Imaging Assembly
    (AIA). SUVI images the Sun in six EUV wavelengths: 9.4, 13.1, 17.1, 19.5,
    28.4, and 30.4 nm. The instrument consists of the main imaging telescope, a
    Guide Telescope (GT) and a Camera Electronics Box (CEB) mechanically integrated
    to the telescope, and a SUVI Electronics Box (SEB). The SEB provides power and
    data interfaces to the spacecraft. The optical chain in the telescope consists
    of thin film entrance filters (metals vapor deposited on a supporting metallic
    mesh), multi-layer coated primary and secondary mirrors, a set of thin film
    analysis filters in two filter wheels near the focal plane, and a charge-coupled
    device (CCD) detector. The CCD consists of 1280 x 1280 pixels with a plate
    scale of 2.5 arcsec and together with the optical system provides a nominal
    53 arcminute square field of view (FOV) from the geostationary orbit. An
    aperture selector with a 60 degree  opening and two multi-segmented mirrors
    enable the imaging in any of the six wavelengths in a single telescope body.
    A nominal 4-minute cadence provides for the observation of the Sun in all
    wavelengths while meeting large dynamic range requirements on single- and
    multi-spectral images. The GT provides Sun-pointing knowledge. In the nominal
    Sun-pointing case, the spacecraft control system uses the GT data to control the
    line-of-sight (LOS) to the Sun.

    Notes
    -----
    SUVI color tables are the same as the AIA color tables for the same
    wavelength with exceptions of SUVI 195 and 284 (which has no direct
    equivalent obviously).  SUVI 195 and 284 images use the AIA 193 & 335 color
    tables respectively.

    References
    ----------
    * `SUVI Instrument Page <https://www.goes-r.gov/spacesegment/suvi.html>`_
    * `Instrument description <https://doi.org/10.3847/2041-8213/aaa28e>`_
    """

    def __init__(self, data, header, **kwargs):

        super().__init__(data, header, **kwargs)

        # Fill in some missing info
        self.meta['detector'] = "SUVI"
        self.meta['telescop'] = "GOES-R"
        self._nickname = self.detector
        self.plot_settings['cmap'] = plt.get_cmap(self._get_cmap_name())
        self.plot_settings['norm'] = ImageNormalize(stretch=source_stretch(self.meta, AsinhStretch(0.01)))

    @property
    def observatory(self):
        """
        Returns the observatory.
        """
        return self.meta['telescop'].split('/')[0]

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume', '').startswith('GOES-R Series Solar Ultraviolet Imager')
