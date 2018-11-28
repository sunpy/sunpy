"""KCor Map subclass definitions"""
# pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from astropy.units import Quantity
from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map import GenericMap
from sunpy.sun import constants
from sunpy.sun import sun
from sunpy.map.sources.source_type import source_stretch
from sunpy.coordinates import get_sunearth_distance

__all__ = ['KCorMap']

class KCorMap(GenericMap):
    """KCor Image Map.

    The COronal Solar Magnetism Observatory (COSMO) K-coronagraph
    (K-Cor) is one of three proposed instruments in the COSMO facility
    suite. It is specifically designed to study the formation and
    dynamics of coronal mass ejections and the evolution of the
    density structure of the low corona. The K-Cor records the
    polarization brightness (pB) formed by Thomson scattering of
    photospheric light by coronal free electrons. The National Center
    for Atmospheric Research (NCAR), via the National Science
    Foundation (NSF), provided full funding for the COSMO K-Cor, which
    was deployed to the Mauna Loa Solar Observatory (MLSO) in Hawaii
    in September 2013, replacing the aging MLSO Mk4 K-coronameter.

    References
    ----------
    * `K-Cor Instrument Page <https://www2.hao.ucar.edu/mlso/instruments/mlso-kcor-coronagraph`_
    """

    def __init__(self, data, header, **kwargs):
        GenericMap.__init__(self, data, header, **kwargs)

        # Fill in some missing info
        self.meta['detector'] = "KCor"
        self._nickname = self.detector
        self.plot_settings['cmap'] = plt.get_cmap(self._get_cmap_name())
        self.plot_settings['norm'] = ImageNormalize(stretch=source_stretch(self.meta, PowerStretch(0.5)))

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an KCor image"""
        return header.get('instrume') == 'EIT'
