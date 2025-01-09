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
    SUIT Image Map

    The Solar Ultraviolet Imaging Telescope (SUIT) is one of the remote sensing
    payloads on board the Aditya-L1 mission of the Indian Space Research Organization (ISRO)
    that was launched on September 02, 2023. SUIT is designed to provide near-simultaneous full-disk
    and region-of-interest images of the Sun at various heights, slicing through the photosphere
    and chromosphere, employing an array of 11 scientifically calibrated filters (3 broad-band & 8 narrow-band)
    strategically positioned within the wavelength range of 200 to 400 nanometers.
    Located at the first Lagrange point, SUIT observes the Sun 24x7, without any interruption.

    The SUIT data is available to download from the PRADAN interface of ISRO.

    References
    ----------
    * `Home Page: <https://suit.iucaa.in/>`_
    * `Mission Page: <https://suit.iucaa.in/about_SUIT>`_
    * `Publications Page: <https://suit.iucaa.in/node/5>`_
    * `Data Download PRADAN Page: <https://pradan.issdc.gov.in/al1>`_
    ----------
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self._nickname = self.detector
        filtername = header.get('FTR_NAME').strip()
        filternorms = {
                "NB01": 0.1,
                "NB02": 9.0,
                "NB03":0.25,
                "NB04":0.2,
                "NB05":2.0,
                "NB06":0.25,
                "NB07":0.25,
                "NB08":0.1,
                "BB01":0.8,
                "BB02":0.4,
                "BB03":0.8,
                }
        self.plot_settings['cmap'] = f"suit_{filtername.lower()}"
        self.plot_settings['title'] = f"SUIT {filtername}:{self.wavelength} - {self.reference_date}"
        self.plot_settings["norm"] = ImageNormalize(stretch=source_stretch(self.meta, AsinhStretch(filternorms.get(filtername, 0.2))), clip=False)

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
        return self.meta.get("MISSION", "Aditya-L1")

    @property
    def instrument(self):
        return "Aditya-L1 (Solar Ultraviolet Imaging Telescope)"

    @property
    def detector(self):
        return self.meta.get("PNAME", "SUIT")

    @property
    def processing_level(self):
        return self.meta.get('F_LEVEL','')

    @property
    def _date_obs(self):
        """
        Time of observation
        """
        return self._get_date('T_OBS')

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
