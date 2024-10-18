from sunpy.net.dataretriever import GenericClient
from scraper import Scraper
__all__ = ["AIASynopticClient"]


class AIASynopticClient(GenericClient):
    """
    A client for retrieving AIA synoptic data from JSOC.

    This client retrieves synoptic AIA data from the following source:
    https://jsoc1.stanford.edu/data/aia/synoptic/

    The synoptic dataset includes lower-resolution (1k) images, with additional
    image processing steps like downsampling and time integration. Characteristics
    of the dataset can be found here: https://jsoc1.stanford.edu/data/aia/synoptic/README.html

    Attributes:
        baseurl (str): The base URL template for retrieving AIA synoptic data.
        pattern (str): The URL pattern for AIA data, including wavelength.
        known_wavelengths (list): A list of known wavelength codes for AIA data.
    """

    known_wavelengths = [171, 193, 211, 304, 335, 1600, 1700]
    url_stem = r"https://jsoc1.stanford.edu/data/aia/synoptic/"
    baseurl = url_stem + r"%Y/%m/%d/H%H00/AIA%Y%m%d_%H%M_(\d{4}).fits"
    pattern = "{}synoptic/{year:4d}/{month:2d}/{day:2d}/H{}/AIA{}_{hour:2d}{minute:2d}_{Wavelength:4d}.fits"

    from sunpy.net import attrs as a
    required = {a.Time, a.Instrument}
    optional = {a.Sample, a.Resolution, a.Wavelength, a.Level}

    @property
    def info_url(self):
        return self.url_stem

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs

        adict = {
            attrs.Instrument: [
                ("AIASynoptic", "AIA Synoptic data from the Atmospheric Imaging Assembly instrument."),
                # ("AIA", "Data from the Atmospheric Imaging Assembly instrument."),
            ],
            attrs.Physobs: [
                ("intensity", "Brightness or intensity of the solar atmosphere at different wavelengths.")
            ],
            attrs.Source: [("SDO", "The Solar Dynamics Observatory.")],
            attrs.Provider: [("JSOC", "Joint Science Operations Center at Stanford.")],
            attrs.Level: [("1.5", "Level 1.5 data processed for specialized analysis.")],
            attrs.Wavelength: [(f"{wv:04d}", f"{wv} AA") for wv in cls.known_wavelengths],
            attrs.Resolution: [("Synoptic", "1024x1024 dataset with 2 minute integration")],
        }
        return adict

    @classmethod
    def pre_search_hook(cls, *args, **kwargs):
        """
        Helper function to return the baseurl, pattern and matchdict
        for the client required by :func:`~sunpy.net.dataretriever.GenericClient.search`
        before using the scraper.
        """
        matchdict = cls._get_match_dict(*args, **kwargs)
        return cls.baseurl, cls.pattern, matchdict
