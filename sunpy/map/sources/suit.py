"""
SUIT Map sub class definitions

"""
from sunpy.map.sources.source_type import source_stretch
import numpy as np
import astropy.units as u
from astropy.coordinates import CartesianRepresentation, HeliocentricMeanEcliptic
from sunpy.map.mapbase import GenericMap
from matplotlib import colors
from astropy.visualization import ImageNormalize
from astropy.visualization import AsinhStretch

__all__ = ['SUITMap']
__author__ = ['Rahul Gopalakrishnan']
__email__ = ['rahulg.astro@gmail.com']

class SUITMap(GenericMap):
    """
    SUIT is blah blah blah

    References:
    -----------
    Blah blah
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self._nickname = self.detector
        self.filter = header.get('FTR_NAME').strip()
        self.wavelnth = str(header.get('WAVELNTH','NA')).strip()
        self.plot_settings['cmap'] = f"suit_{self.filter.lower()}"
        self.plot_settings['title'] = f"SUIT {self.filter}:{self.wavelnth} A - {self.reference_date}"
        norm = self.suit_norm()

        if self.filter in norm:
            self.plot_settings["norm"] = norm[self.filter]
        else:
            self.plot_settings["norm"] = ImageNormalize(
                stretch=source_stretch(self.meta, AsinhStretch(0.01)), clip=False
            )

    def suit_norm(self):
        a = np.array(self.data)
        sliced_array = a[(a > 2000) & (a < 50000)]
        baseline = np.median(sliced_array) if len(sliced_array) > 0 else 1000
        norm_params = {
            'NB05': (0.8, 0.2 * baseline, 3 * baseline),
            'NB03': (0.8, 0.2 * baseline, 3 * baseline),
            'NB04': (0.9, 0.2 * baseline, 2.5 * baseline),
            'NB02': (0.9, 0.2 * baseline, 3 * baseline),
            'NB07': (0.8, 0.2 * baseline, 2.5 * baseline),
            'NB06': (0.8, 0.2 * baseline, 2.5 * baseline),
            'BB03': (0.8, 0.2 * baseline, 2.5 * baseline),
            'BB02': (0.8, 0.2 * baseline, 2.3 * baseline),
            'BB01': (.8, 0.2 * baseline, 2.5 * baseline),
            'NB01': (0.9, 0.1 * baseline, 3 * baseline),
            'NB08': (0.8, 0.2 * baseline, 3 * baseline),
        }

        norm = {key: colors.PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax)
                 for key, (gamma, vmin, vmax) in norm_params.items()}
        return norm

    @property
    def _supported_observer_coordinates(self):
        return [(('haex_obs', 'haey_obs', 'haez_obs'), {'x': self.meta.get('haex_obs'),
                                                        'y': self.meta.get('haey_obs'),
                                                        'z': self.meta.get('haez_obs'),
                                                        'unit': u.m,
                                                        'representation_type': CartesianRepresentation,
                                                        'frame': HeliocentricMeanEcliptic})
                ] + super()._supported_observer_coordinates

    @property
    def observatory(self):
        return 'Aditya-L1'

    @property
    def instrument(self):
        return 'Aditya-L1 (Solar Ultraviolet Imaging Telescope)'

    @property
    def detector(self):
        return 'SUIT'

    @property
    def processing_level(self):
        return self.meta.get('F_LEVEL','')

    @property
    def date_obs(self):
        return self.meta.get('T_OBS','')

    @property
    def reference_date(self):
        """
        The reference time is the time of Shutter open time. Does not include the exposure time.
        """
        return self._get_date('T_OBS') or super().reference_date


    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if header corresponds to a SUIT image

        """
        payload = header.get('PNAME', '').strip().lower() == 'suit'
        return  payload
