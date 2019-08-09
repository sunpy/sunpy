"""SUVI Map subclass definitions"""
from __future__ import absolute_import, print_function, division

# pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

import matplotlib.pyplot as plt

from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import AsinhStretch

from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ["SUVIMap"]


class SUVIMap(GenericMap):
    """SUVI Image Map.

    The Solar Ultraviolet Imager (SUVI) is a normal-incidence Cassegrain EUV
    telescope onboard the latest of the Geostationary Operational Environmental
    Satellite (GOES) missions (GOES-16, formerly known as GOES-R).
    It is similar to Atmospheric Imaging Assembly (AIA). It operates in
    geostationary orbit above the Americas at 75.2 degree W. It's primary
    purpose is to support NOAA's goal to characterize solar features and detect
    events that lead to space weather. It uses a filter wheel to image the Sun
    in six EUV wavelength corresponding to known coronal emission lines:

    - 9.4 nm (FeXVIII)
    - 13.1 nm (FeXXI)
    - 17.1 nm (FeIX/X)
    - 19.5 nm (FeXII)
    - 28.4 nm (FeXV)
    - 30.4 nm (HeII)

    The focal plane consists of a CCD detector with 1280 x 1280 pixels. The
    plate scale is 2.5 arcsec per pixel. The field of view is therefore almost
    twice the size of the Sun (53 arcmin) and extends out to 1.6 solar radii in
    the horizontal direction and 2.3 solar radii in the diagonal. It provides
    observations in each wavelength at multiple exposure times every 4 minutes.

    It began operating on ???.

    Notes
    -----
    SUVI uses the same color tables as AIA for the matching wavelengths.
    SUVI 195 and 284 images use the AIA 193 & 335 color tables respectively.

    References
    ----------
    * `GOES-R Mission<https://www.goes-r.gov>`_
    * `SUVI Instrument Page <https://www.goes-r.gov/spacesegment/suvi.html>`_
    * `GOES-16 on wikipedia <https://en.wikipedia.org/wiki/GOES-16>`_
    * `Instrument paper: Coronal Imaging with the Solar UltraViolet Image <http://doi.org/10.1007/s11207-019-1411-0>`_
    * `User's Guide <https://www.goes-r.gov/users/docs/PUG-L1b-vol3.pdf>`_
    * `Level 1b Readme <https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l1b/suvi-l1b-fe094/ReadMe.pdf>`_
    * `Data archive <https://www.ngdc.noaa.gov/stp/satellite/goes-r.html>`_,
      `Level 1b data <https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l1b/>`_,
      `Level 2 data <https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/>`_
    """

    def __init__(self, data, header, **kwargs):

        super().__init__(data, header, **kwargs)

        # Fill in some missing info
        self.meta["detector"] = "SUVI"
        self.meta["telescop"] = "GOES-R"
        self._nickname = self.detector
        self.plot_settings["cmap"] = plt.get_cmap(self._get_cmap_name())
        self.plot_settings["norm"] = ImageNormalize(
            stretch=source_stretch(self.meta, AsinhStretch(0.01))
        )

    @property
    def observatory(self):
        """
        Returns the observatory.
        """
        return self.meta["telescop"].split("/")[0]

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an AIA image"""
        return header.get("instrume", "").startswith(
            "GOES-R Series Solar Ultraviolet Imager"
        )
