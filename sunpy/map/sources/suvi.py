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

    The GOES-R Series Solar Ultraviolet Imager (SUVI)....(description of the
    instrument).

    SUVI color tables are the same as the AIA color tables for the same
    wavelength with exceptions of SUVI 195 and 284 (which has no direct
    equivalent obviously).  SUVI 195 and 284 images use the AIA 193 & 335 color
    tables respectively.

    References
    ----------
    * `SUVI Mission Page <https://sdo.gsfc.nasa.gov/>`_
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

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
