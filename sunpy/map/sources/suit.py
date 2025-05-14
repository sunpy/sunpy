"""
SUIT map subclass definitions
"""
from astropy.visualization import AsinhStretch, ImageNormalize

from sunpy.map.mapbase import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ["SUITMap"]

class SUITMap(GenericMap):
    """
    SUIT Image Map.

    The Solar Ultraviolet Imaging Telescope (SUIT) is one of the remote sensing
    payloads on board the Aditya-L1 mission of the Indian Space Research Organization (ISRO)
    that was launched on September 2, 2023 and located at the first Lagrange point.
    SUIT is designed to provide near-simultaneous full-disk
    and region-of-interest images of the Sun at various heights, slicing through the photosphere
    and chromosphere, employing an array of 11 scientifically calibrated filters (3 broad-band and 8 narrow-band)
    strategically positioned within the wavelength range of 200 to 400 nanometers.

    The SUIT data is available to download from the PRADAN interface of ISRO.

    References
    ----------
    * `Home Page: <https://suit.iucaa.in/>`__
    * `Mission Page: <https://suit.iucaa.in/about_SUIT>`__
    * `Publications Page: <https://suit.iucaa.in/node/5>`__
    * `Data Download PRADAN Page: <https://pradan.issdc.gov.in/al1>`__
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self._nickname = self.detector
        filtername = header.get("FTR_NAME", "").strip()
        filternorms = {
            "NB01": 0.1,
            "NB02": 9.0,
            "NB03": 0.25,
            "NB04": 0.2,
            "NB05": 2.0,
            "NB06": 0.25,
            "NB07": 0.25,
            "NB08": 0.1,
            "BB01": 0.8,
            "BB02": 0.4,
            "BB03": 0.8,
        }
        self.plot_settings["cmap"] = f"suit_{filtername.lower()}"
        self.plot_settings["norm"] = ImageNormalize(stretch=source_stretch(self.meta, AsinhStretch(filternorms.get(filtername, 0.2))), clip=False)

    @property
    def observatory(self):
        return self.meta.get("MISSION", "Aditya-L1")

    @property
    def instrument(self):
        return "SUIT"

    @property
    def detector(self):
        return self.meta.get("PNAME", "SUIT")

    @property
    def processing_level(self):
        return self.meta.get("F_LEVEL", "")

    @property
    def _date_obs(self):
        """
        Time of observation
        """
        return self._get_date("T_OBS")

    @property
    def reference_date(self):
        """
        The reference time is the time of Shutter open time. Does not include the exposure time.
        """
        return self._get_date("T_OBS") or super().reference_date

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if header corresponds to a SUIT image
        """
        return header.get("PNAME", "").strip().lower() == "suit"
