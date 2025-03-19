"""STEREO Map subclass definitions"""

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"


import astropy.units as u
from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ['EUVIMap', 'CORMap', 'HIMap']


class EUVIMap(GenericMap):
    """STEREO-SECCHI EUVI Image Map

    EUVI is an extreme ultraviolet (EUV) imager. Part of the STEREO-SECCHI
    suite it observes the Sun from 1 to 1.7 solar radii. It is capable of
    observing at 304 (He II), 171 (Fe IX), 195 (Fe XII), and 284 (Fe XV)
    Angstroms.

    References
    ----------
    * `STEREO Mission Page <https://stereo.gsfc.nasa.gov/>`__
    * `STEREO SECCHI <https://secchi.nrl.navy.mil/>`__
    * `Instrument Page <http://secchi.lmsal.com/EUVI/>`__
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self._nickname = f"{self.detector}-{self.observatory[-1]}"
        self.plot_settings['cmap'] = f'euvi{int(self.wavelength.value):d}'
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, PowerStretch(0.25)), clip=False)

    def _rotation_matrix_from_crota(self):
        return super()._rotation_matrix_from_crota('CROTA')

    @property
    def waveunit(self):
        unit = self.meta.get("waveunit", "Angstrom")
        return u.Unit(unit)

    @property
    def rsun_arcseconds(self):
        return self.meta.get('rsun', None)

    @property
    def rsun_obs(self):
        rsun_arcseconds = self.rsun_arcseconds

        if rsun_arcseconds is None:
            rsun_arcseconds = super().rsun_obs

        return u.Quantity(rsun_arcseconds, 'arcsec')

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EUVI image"""
        return header.get('detector') == 'EUVI'


class CORMap(GenericMap):
    """STEREO-SECCHI CORonograph Image Map.

    Part of the STEREO-SECCHI suite of remote sensing telescopes,
    COR is a set of two coronographs (COR1, COR2) onboard STEREO.
    They are both traditional Lyot coronagraphs.

    The COR1 detectors observes from 1.3 to 4 solar radii while the
    COR2 detectors observe a range from 2 to 15 solar radii.

    References
    ----------
    * `STEREO Mission Page <https://stereo.gsfc.nasa.gov/>`__
    * `STEREO SECCHI <https://secchi.nrl.navy.mil/>`__
    * `COR1 Instrument Page <https://cor1.gsfc.nasa.gov>`__
    * `COR2 Instrument Page <https://secchi.nrl.navy.mil//index.php?p=cor2>`__
    * `COR1 User Guide <https://cor1.gsfc.nasa.gov/guide/>`__
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self._nickname = f"{self.detector}-{self.observatory[-1]}"
        self.plot_settings['cmap'] = f'stereocor{self.detector[-1]!s}'
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, PowerStretch(0.5)), clip=False)

    @property
    def measurement(self):
        # TODO: This needs to do more than white-light. Should give B, pB, etc.
        return "white-light"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an COR image"""
        return str(header.get('detector', '')).startswith('COR')


class HIMap(GenericMap):
    """STEREO-SECCHI Heliospheric Imager (HI) Map.

    The HI is a wide-angle visible-light imaging system
    for the detection of coronal mass ejection (CME) events
    in interplanetary space and, in particular, of events
    directed towards the Earth.

    The Heliospheric imager consists of two instruments, the HI-1 and HI-2.
    The HI1 observes from 15-80 solar radii while HI2 observes from 80-215
    solar radii.

    References
    ----------
    * `STEREO Mission Page <https://stereo.gsfc.nasa.gov/>`__
    * `STEREO SECCHI <https://secchi.nrl.navy.mil>`__
    * `HI Instrument Page <http://www.stereo.rl.ac.uk>`__
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self._nickname = f"{self.detector}-{self.observatory[-1]}"
        self.plot_settings['cmap'] = f'stereohi{self.detector[-1]!s}'
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, PowerStretch(0.25)), clip=False)

    @property
    def measurement(self):
        # TODO: This needs to do more than white-light. Should give B, pB, etc.
        return "white-light"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an COR image"""
        return str(header.get('detector', '')).startswith('HI')
