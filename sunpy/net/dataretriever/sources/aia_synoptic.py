import logging
from sunpy.net import attrs, Fido
from sunpy.net.attr import SimpleAttr
from sunpy.net.dataretriever import GenericClient
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.attrs.aia_synoptic import AIASynopticData
import astropy.units as u
import os
from datetime import timedelta, datetime


# Logger setup
def setup_logger():
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

    baseurl = r"https://jsoc1.stanford.edu/data/aia/synoptic/%Y/%m/%d/H%H00/"
    pattern = "{baseurl}AIA%Y%m%d_%H%M_{wavelength}.fits"
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
    def _can_handle_query(cls, *query: SimpleAttr) -> bool:
        """
        Determines if the client can handle the given query.

        Parameters:
            query (SimpleAttr): A set of query attributes.

        Returns:
            bool: True if the client can handle the query, False otherwise.
        """
        required_attrs = {attrs.Time, AIASynopticData}
        all_attrs = {type(x) for x in query}

        return required_attrs.issubset(all_attrs) and any(isinstance(x, AIASynopticData) for x in query)

    def search(self, *query: SimpleAttr) -> QueryResponse:
        """
        Perform a search query for AIA synoptic data.

        Parameters:
            query (SimpleAttr): The query parameters including time, wavelength, and cadence.

        Notes:
            If AIASynopticData is present, resolution defaults to 1k.
            If a resolution is specified alongside AIASynopticData, it will be overridden to 1k.

        Returns:
            QueryResponse: A response object containing the result of the query.
        """
        time_range = None
        wavelengths = []
        cadence_seconds = None
        use_synoptic_data = False

        for q in query:
            if isinstance(q, attrs.Time):
                time_range = q
            elif isinstance(q, attrs.Wavelength):
                wavelengths.append(int(q.min.to(u.angstrom).value))
            elif isinstance(q, attrs.Sample):
                cadence_seconds = q.value
            elif isinstance(q, AIASynopticData):
                use_synoptic_data = True
                # If synoptic data is used, enforce 1k resolution
                if any(isinstance(attr, attrs.Resolution) for attr in query):
                    logger.warning("Resolution is overridden to 1k due to the use of AIASynopticData.")

        if not time_range:
            logger.error("Time range must be specified for the AIASynopticClient.")
            raise ValueError("Time range must be specified for the AIASynopticClient.")

        wavelengths = [str(wl).zfill(4) for wl in (wavelengths or self.known_wavelengths)]
        urls = self._generate_urls(time_range, wavelengths, cadence_seconds)
        return self._prepare_query_response(urls)

    def _prepare_query_response(self, urls: list) -> QueryResponse:
        """
        Prepare the query response by populating the necessary fields.

        Parameters:
            urls (list): A list of URLs corresponding to the requested data.

        Returns:
            QueryResponse: A response object containing the prepared query response.
        """
        from sunpy.net.dataretriever.client import QueryResponseTable

        data = {
            "Start Time": [],
            "Instrument": ["AIA"] * len(urls),
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

    def _generate_urls(self, time_range: attrs.Time, wavelengths: list, cadence_seconds: int = None) -> list:
        """
        Generate a list of URLs for AIA synoptic data given a time range, wavelengths, and cadence.

        Parameters:
            time_range (attrs.Time): The time range for the query.
            wavelengths (list): List of wavelength values.
            cadence_seconds (int, optional): The cadence in seconds between each time step.

        Returns:
            list: URLs corresponding to the requested data.
        """
        urls = []
        current_time = time_range.start.datetime
        end_time = time_range.end.datetime
        cadence = timedelta(seconds=cadence_seconds) if cadence_seconds else timedelta(minutes=1)

        logger.info(f"Using cadence: {cadence}.")  # Log cadence

        while current_time <= end_time:
            baseurl = current_time.strftime(self.baseurl)
            for wavelength in wavelengths:
                urls.append(self.pattern.format(baseurl=baseurl, wavelength=str(wavelength).zfill(4)))
            current_time += cadence

        logger.debug(f"Generated {len(urls)} URLs for download.")
        return urls

    def fetch(self, query_result, *, path: str, downloader, **kwargs):
        """
        Fetch the data for a given query result and download it to the specified path.

        Parameters:
            query_result: The result from the query to fetch.
            path (str): The path where the data will be downloaded.
            downloader: The downloader object responsible for fetching files.
            kwargs: Additional keyword arguments for configuration.
        """
        try:
            download_path = os.path.dirname(path)
            downloader.max_conn = kwargs.get("max_conn", 10)
            for record in query_result:
                downloader.enqueue_file(record["url"], path=download_path)
            return downloader.download()
        except (IOError, OSError) as e:
            logger.error(f"File error while fetching data: {e}")
            raise
        except Exception as e:
            logger.error(f"Error while fetching data: {e}")
            raise


# Register the client with Fido
if AIASynopticClient not in Fido.registry:
    Fido.registry[AIASynopticClient] = AIASynopticClient._can_handle_query
    logger.info("Synoptic Fido Client Loaded!")
