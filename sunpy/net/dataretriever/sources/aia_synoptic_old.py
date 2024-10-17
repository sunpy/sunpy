import os
import logging
from datetime import datetime, timedelta
import astropy.units as u
from sunpy.net import attrs
from sunpy.net.dataretriever import GenericClient
from sunpy.net.dataretriever.attrs.aia_synoptic import AIASynopticData
from sunpy.net.dataretriever.client import QueryResponse


# Logger setup
def setup_logger() -> logging.Logger:
    logger = logging.getLogger(__name__)
    log_level = os.getenv("LOG_LEVEL", "INFO").upper()
    logger.setLevel(log_level)

    if not logger.hasHandlers():
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger


logger = setup_logger()

__all__ = ["AIASynopticClient"]


class AIASynopticClient(GenericClient):
    """
    A client for retrieving AIA synoptic data from JSOC.

    This client retrieves synoptic AIA data from the following source:
    https://jsoc1.stanford.edu/data/aia/synoptic/.

    The synoptic dataset consists of lower-resolution (1k) images, with level 1.5
    image processing steps like downsampling and time integration. Characteristics
    of the dataset can be found here: https://jsoc1.stanford.edu/data/aia/synoptic/README.html.

    Attributes:
        baseurl (str): The base URL template for retrieving AIA synoptic data.
        pattern (str): The URL pattern for AIA data, including wavelength.
        known_wavelengths (list): A list of known wavelength codes for AIA data.
    """

    baseurl = r"https://jsoc1.stanford.edu/data/aia/synoptic/%Y/%m/%d/H%H00/"
                r""
    pattern = "{}AIA%Y%m%d_%H%M_{wavelength}.fits"
    known_wavelengths = [
        "0094",
        "0131",
        "0171",
        "0193",
        "0211",
        "0304",
        "1600",
        "1700",
    ]

    @classmethod
    def _attrs_module(cls):
        return "aia_synoptic", "sunpy.net.dataretriever.attrs.aia_synoptic"

    @classmethod
    def _can_handle_query(cls, *query) -> bool:
        """Determines if the client can handle the given query."""
        synoptic = attrs.Instrument("AIASynoptic")
        required_attrs = {attrs.Time, (synoptic | AIASynopticData)}
        all_attrs = {type(attr) for attr in query}
        return required_attrs.issubset(all_attrs) and (
            any(isinstance(attr, AIASynopticData) for attr in query)
            or any(isinstance(attr, synoptic) for attr in query)
        )

    def search(self, *query) -> QueryResponse:
        """Perform a search query for AIA synoptic data."""
        time_range = None
        wavelengths = []
        cadence_seconds = None

        for q in query:
            if isinstance(q, attrs.Time):
                time_range = q
            elif isinstance(q, attrs.Wavelength):
                wl_value = int(q.min.to(u.angstrom).value)
                if str(wl_value).zfill(4) in self.known_wavelengths:
                    wavelengths.append(wl_value)
                else:
                    logger.warning(f"Unknown wavelength {wl_value} ignored.")
            elif isinstance(q, attrs.Sample):
                cadence_seconds = q.value
            elif isinstance(q, AIASynopticData):
                logger.warning("Resolution is overridden to 1k due to the use of AIASynopticData.")

        if not time_range:
            logger.error("Time range must be specified for the AIASynopticClient.")
            raise ValueError("Time range must be specified for the AIASynopticClient.")

        wavelengths = [str(wl).zfill(4) for wl in (wavelengths or self.known_wavelengths)]
        urls = self._generate_urls(time_range, wavelengths, cadence_seconds)
        return self._prepare_query_response(urls)

    def _generate_urls(self, time_range: attrs.Time, wavelengths: list, cadence_seconds: int = None) -> list:
        """Generate a list of URLs for AIA synoptic data."""
        urls = []
        current_time = time_range.start.datetime
        end_time = time_range.end.datetime
        cadence = timedelta(seconds=cadence_seconds) if cadence_seconds else timedelta(minutes=1)

        while current_time <= end_time:
            baseurl = current_time.strftime(self.baseurl)
            pattern = current_time.strftime(self.pattern)
            for wavelength in wavelengths:
                urls.append(pattern.format(baseurl=baseurl, wavelength=str(wavelength).zfill(4)))
            current_time += cadence

        logger.debug(f"Generated {len(urls)} URLs for download.")
        return urls

    def _prepare_query_response(self, urls: list) -> QueryResponse:
        """Prepare the query response by populating the necessary fields."""
        from sunpy.net.dataretriever.client import QueryResponseTable

        data = {
            "Start Time": [],
            "Instrument": ["AIASynoptic"] * len(urls),
            # "Provider": ["JSOCweb"] * len(urls),
            "Resolution": ["1024"] * len(urls),
            "Wavelength": [],
            "url": [],
            "fileid": [],
        }

        for url in urls:
            filename = os.path.basename(url)
            parts = filename[3:-5].split("_")
            if len(parts) == 3:
                try:
                    start_time = datetime.strptime(parts[0] + parts[1], "%Y%m%d%H%M")
                except ValueError:
                    start_time = None
                data["Start Time"].append(start_time)
                data["Wavelength"].append(int(parts[2]))
            else:
                data["Start Time"].append(None)
                data["Wavelength"].append(None)
            data["url"].append(url)
            data["fileid"].append(filename)

        return QueryResponse(QueryResponseTable(data, client=self))


# Register the client with Fido


def main():
    start_time = "2020-01-01T00:00:00.00"
    end_time = "2020-01-01T23:59:59.999"

    # Client = AIASynopticClient()
    from sunpy.net import Fido, attrs as a

    # from sunpy.net import Fido
    # from sunpy.net.dataretriever.sources.aia_synoptic import AIASynopticClient

    # Register the AIASynopticClient with Fido
    # Fido.register_client(AIASynopticClient)

    request = Fido.search(
        attrs.Instrument("AIA"),
        attrs.Time(start_time, end_time),
        attrs.Sample(1 * u.hour),
        attrs.Wavelength(171 * u.angstrom),
        attrs.ExtentType("SYNOPTIC"),
    )
    print(request)


if __name__ == "__main__":
    main()
