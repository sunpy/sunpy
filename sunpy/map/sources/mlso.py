import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ['KCorMap']


class KCorMap(GenericMap):
    """
    K-Cor Image Map.

    The COronal Solar Magnetism Observatory (COSMO) K-coronagraph (K-Cor) is one of three proposed
    instruments in the COSMO facility suite. It is specifically designed to study the formation
    and dynamics of coronal mass ejections and the evolution of the density structure of
    the low corona. The K-Cor records the polarization brightness (pB) formed by Thomson scattering
    of photospheric light by coronal free electrons. The National Center for Atmospheric
    Research (NCAR), via the National Science Foundation (NSF), provided full funding for the
    COSMO K-Cor, which was deployed to the Mauna Loa Solar Observatory (MLSO) in Hawaii in
    September 2013, replacing the aging MLSO Mk4 K-coronameter.

    Notes
    -----
    Observer location: The standard K-Cor metadata does not include the full 3D
    location of the observer. There are 2D Carrington heliographic coordinates, but
    using them is not recommended because the calculation of Carrington longitude
    differs from ``sunpy`` (see :ref:`sunpy-topic-guide-coordinates-carrington`).
    Instead, we assume the default observer location to be the geographic location of MLSO.

    References
    ----------
    * `COSMO Mission Page <https://www2.hao.ucar.edu/cosmo>`__
    * `KCOR Instrument Page <https://www2.hao.ucar.edu/mlso/instruments/mlso-kcor-coronagraph>`__
    """

    # MLSO location per Wikipedia (https://en.wikipedia.org/wiki/Mauna_Loa_Solar_Observatory)
    _earth_location = EarthLocation(-155.576*u.deg, 19.536*u.deg, 3394*u.m)

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self._nickname = self.detector

        self.plot_settings['cmap'] = self._get_cmap_name()
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, PowerStretch(0.25)), clip=False)
        # Negative value pixels can appear that lead to ugly looking images.
        # This can be fixed by setting the lower limit of the normalization.
        self.plot_settings['norm'].vmin = 0.0

    def _get_cmap_name(self):
        """Build the default color map name."""
        return self.detector.lower()

    @property
    def observatory(self):
        return "MLSO"

    @property
    def detector(self):
        return "KCor"

    @property
    def waveunit(self):
        """
        If the WAVEUNIT FITS keyword is not present, defaults to nanometers.
        """
        unit = self.meta.get("waveunit", "nm")
        return u.Unit(unit)

    @property
    def _default_observer_coordinate(self):
        return SkyCoord(self._earth_location.get_itrs(self.date)).heliographic_stonyhurst

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to a COSMO image"""
        return header.get('instrume') == 'COSMO K-Coronagraph'
