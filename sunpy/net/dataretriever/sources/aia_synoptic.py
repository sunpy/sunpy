from sunpy.net.dataretriever import GenericClient
from sunpy.net.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ["AIASynopticClient"]


class AIASynopticClient(GenericClient):
    """
    A client for retrieving AIA synoptic data from JSOC.

    This client retrieves synoptic AIA data from the following source:
    https://jsoc1.stanford.edu/data/aia/synoptic/

    The synoptic dataset includes lower-resolution (1k) images, with additional
    image processing steps like downsampling and time integration. Characteristics
    of the dataset can be found here: https://jsoc1.stanford.edu/data/aia/synoptic/README.html

    - If AIASynopticData is present, resolution defaults to 1k and any user-specified
      resolution will be overridden.

    Attributes:
        baseurl (str): The base URL template for retrieving AIA synoptic data.
        pattern (str): The URL pattern for AIA data, including wavelength.
        known_wavelengths (list): A list of known wavelength codes for AIA data.
    """

    # example_url = info_url + r"2022/12/06/H0000/AIA20221206_0000_0171.fits"
    known_wavelengths = [
        171,
        193,
        211,
        304,
        335,
        1600,
        1700,
    ]
    info_url = r"https://jsoc1.stanford.edu/data/aia/synoptic/"
    baseurl = info_url + r"%Y/%m/%d/H%H00/AIA%Y%m%d_%H%M_.....fits"
    pattern = info_url + r"%Y/%m/%d/H%H00/AIA%Y%m%d_%H%M_.....fits"
    patty = "{}synoptic/{year:4d}/{month:2d}/{day:2d}/H{}/AIA{}_{hour:2d}{minute:2d}_{Wavelength:4d}.fits"
    # pattern = r"https://jsoc1.stanford.edu/data/aia/synoptic/(\d{4})/(\d{2})/(\d{2})/H(\d{2})00/AIA(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})_(?P<hour>\d{2})(?P<minute>\d{2})_(?P<wavelength>\d{4})\.fits"

    # pattern = extractor = (
    #     "https://jsoc1.stanford.edu/data/aia/synoptic/{year}/{month}/{day}/H{hour}00/AIA{year_full}{month_full}{day_full}_{hour_full}{minute}_{wavelength}.fits"
    # )

    from sunpy.net import attrs as a

    required = {a.Time, a.Instrument}

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs

        adict = {
            attrs.Instrument: [
                ("AIASynoptic", "AIA Synoptic data from the Atmospheric Imaging Assembly instrument.")
            ],
            attrs.Physobs: [
                ("intensity", "Brightness or intensity of the solar atmosphere at different wavelengths.")
            ],
            attrs.Source: [("SDO", "The Solar Dynamics Observatory.")],
            attrs.Provider: [("JSOC", "Joint Science Operations Center at Stanford.")],
            attrs.Level: [("1.5", "Level 1.5 data processed for specialized analysis.")],
            attrs.Wavelength: [(f"{wv:04d}", f"{wv} AA") for wv in cls.known_wavelengths],
        }
        return adict


# Register values upon module import
if __name__ == "__main__":
    from sunpy.net import Fido
    from sunpy.net import attrs as a
    import astropy.units as u

    AIASynopticClient.register_values()  # Ensure registration is done

    # Example usage
    time_range = a.Time("2023-10-11 00:00:00", "2023-10-11 06:00:00")
    instrument = a.Instrument("AIASynoptic")
    # wavelength = a.Wavelength(193 * u.angstrom)
    # sample = a.Sample(4 * u.minute)
    results = Fido.search(time_range, instrument)  # , wavelength)

    print(results)

    # @classmethod
    # def pre_search_hook(cls, *args, **kwargs):
    #     """
    #     Helper function to return the baseurl, pattern and matchdict
    #     for the client required by :func:`~sunpy.net.dataretriever.GenericClient.search`
    #     before using the scraper.
    #     """
    #     matchdict = cls._get_match_dict(*args, **kwargs)
    #     return cls.baseurl, cls.pattern, matchdict

    # def _can_handle_query(self, *query):
    #     from sunpy.net import attrs as a

    #     required = {a.Instrument, a.Wavelength}
    #     all_attrs = {type(x) for x in query}
    #     if not required.issubset(all_attrs):
    #         return False
    #     for x in query:
    #         if isinstance(x, a.Instrument) and x.value.lower() != "aiasynoptic":
    #             return False
    #         # if isinstance(x, a.Wavelength) and x.unit != "angstrom":
    #         #     return False
    #     return True
