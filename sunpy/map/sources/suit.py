"""
SUIT Map sub class definitions

"""
from sunpy.map.sources.source_type import source_stretch
import numpy as np
import astropy.units as u
from astropy.coordinates import CartesianRepresentation, HeliocentricMeanEcliptic
from sunpy.map.mapbase import GenericMap
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
        self.plot_settings['cmap'] = f"suit_{self.filter.lower()}"
        self.plot_settings['title'] = f"SUIT {self.filter} - {self.reference_date}"
        self.plot_settings["norm"] = ImageNormalize(
            stretch=source_stretch(self.meta, AsinhStretch(0.01)), clip=False
        )

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
