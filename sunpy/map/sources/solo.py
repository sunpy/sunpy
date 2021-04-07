"""Solar Orbiter Map subclass definitions"""
import astropy.units as u
from astropy.coordinates import CartesianRepresentation
from astropy.visualization import AsinhStretch, ImageNormalize

from sunpy.coordinates import HeliocentricInertial
from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ['EUIMap']


class EUIMap(GenericMap):
    """
    EUI Image Map

    The Extreme Ultraviolet Imager (EUI) is a remote sensing instrument onboard the
    Solar Orbiter (SolO) spacecraft. EUI has three telescopes that image the Sun in
    Lyman-alpha (1216 Å) and the EUV (174 Å and 304 Å). The three telescopes are the
    Full Sun Imager (FSI) and two High Resolution Imagers (HRI). The FSI images the
    whole Sun in both 174 Å and 304 Å. The EUV and Lyman-alpha HRI telescopes image a
    1000"-by-1000" patch in 174 Å and 1216 Å, respectively.

    References
    ----------
    * `Solar Orbiter Mission Page <https://sci.esa.int/web/solar-orbiter/>`_
    * `EUI Instrument Page <https://wwwbis.sidc.be/EUI/EUI/EUI/EUI/EUI/>`_
    * `Instrument Paper <https://doi.org/10.1051/0004-6361/201936663>`_
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self._nickname = self.detector
        self.plot_settings['cmap'] = self._get_cmap_name()
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, AsinhStretch(0.01)), clip=False)
        self.meta['exptime'] = self.meta.get('xposure', 0.0)
        # The level number is prepended by the letter L
        if self.meta.get('level'):
            self.meta['lvl_num'] = int(self.meta.get('level')[1:])
        # DN is not a FITS standard unit, so convert to counts
        if self.meta.get('bunit', None) == 'DN':
            self.meta['bunit'] = 'ct'
        if self.meta.get('bunit', None) == 'DN/s':
            self.meta['bunit'] = 'ct/s'

    @property
    def _supported_observer_coordinates(self):
        return [(('hcix_obs', 'hciy_obs', 'hciz_obs'),
                 {'x': self.meta.get('hcix_obs'),
                  'y': self.meta.get('hciy_obs'),
                  'z': self.meta.get('hciz_obs'),
                  'unit': u.m,
                  'representation_type': CartesianRepresentation,
                  'frame': HeliocentricInertial})] + super()._supported_observer_coordinates

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EUI image"""
        is_solo = 'solar orbiter' in str(header.get('obsrvtry', '')).lower()
        is_eui = str(header.get('instrume', '')).startswith('EUI')
        return is_solo and is_eui
