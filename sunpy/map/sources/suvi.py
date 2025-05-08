"""SUVI Map subclass definitions"""
import astropy.units as u
from astropy.coordinates import CartesianRepresentation
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"


__all__ = ["SUVIMap"]


class SUVIMap(GenericMap):
    """
    SUVI Image Map.

    The Solar Ultraviolet Imager (SUVI) is a normal-incidence Cassegrain EUV
    telescope on board the latest of the Geostationary Operational Environmental
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

    GOES-16 was launched on November 16, 2016, and became operational as NOAA's
    GOES East on December 18, 2017, replacing GOES-13.

    Notes
    -----
    SUVI uses the same color tables as AIA for the matching wavelengths.
    SUVI 195 and 284 images use the AIA 193 & 335 color tables respectively.

    Observer location: We use the ECEF coordinates provided in the FITS header for the spacecraft
    location even when coordinates in other frames are provided due to accuracy concerns over the
    coordinate transformations used in the SUVI data pipeline. There could still be a small
    discrepancy because the definition of the ECEF frame used by SUVI may not exactly match the
    definition of the ITRS frame used by SunPy to interpret the header values.

    Note that some Level 1b files cannot be loaded due to errors in the header.

    References
    ----------
    * `GOES-R Mission <https://www.goes-r.gov>`__
    * `SUVI Instrument Page <https://www.goes-r.gov/spacesegment/suvi.html>`__
    * `GOES-16 on Wikipedia <https://en.wikipedia.org/wiki/GOES-16>`__
    * Recommended instrument paper: :cite:t:`seaton_observations_2018`
    * `User's Guide <https://www.goes-r.gov/users/docs/PUG-L1b-vol3.pdf>`__
    * `Level 1b Readme <https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l1b/suvi-l1b-fe094/ReadMe.pdf>`__
    * `Data archive <https://www.ngdc.noaa.gov/stp/satellite/goes-r.html>`__
    * `Level 1b data <https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l1b/>`__
    * `Level 2 data <https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/>`__
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self._nickname = self.detector
        self.plot_settings["cmap"] = self._get_cmap_name()
        self.plot_settings["norm"] = ImageNormalize(
            stretch=source_stretch(self.meta, AsinhStretch(0.01)), clip=False
        )

    @property
    def _supported_observer_coordinates(self):
        return [(('obsgeo-x', 'obsgeo-y', 'obsgeo-z'), {'x': self.meta.get('obsgeo-x'),
                                                        'y': self.meta.get('obsgeo-y'),
                                                        'z': self.meta.get('obsgeo-z'),
                                                        'unit': u.m,
                                                        'representation_type': CartesianRepresentation,
                                                        'frame': "itrs"})
                ] + super()._supported_observer_coordinates

    @property
    def observatory(self):
        return "GOES-R"

    @property
    def detector(self):
        return "SUVI"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an AIA image"""
        return str(header.get("instrume", "")).startswith(
            "GOES-R Series Solar Ultraviolet Imager"
        )
