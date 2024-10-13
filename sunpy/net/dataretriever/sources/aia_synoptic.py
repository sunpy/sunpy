import logging
from sunpy.net import attrs, Fido
from sunpy.net.attr import SimpleAttr
from sunpy.net.dataretriever import GenericClient
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.attrs.aia_synoptic import AIASynopticData
import astropy.units as u
import os
from datetime import timedelta, datetime

# Setup the logger
logger = logging.getLogger(__name__)
logger.setLevel(os.getenv("LOG_LEVEL", "INFO"))  # Allow log level to be configurable

# Create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# Create formatter and add it to the handler
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)

# Add the handler to the logger
logger.addHandler(ch)

__all__ = ["AIASynopticClient"]


class AIASynopticClient(GenericClient):
    """
    A client for retrieving AIA synoptic data from JSOC.

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
        optional_attrs = {attrs.Instrument, attrs.Wavelength, attrs.Sample}
        all_attrs = {type(x) for x in query}

        if not required_attrs.issubset(all_attrs):
            return False

        has_synoptic_data_attr = any(isinstance(x, AIASynopticData) for x in query)
        return has_synoptic_data_attr

    def search(self, *query: SimpleAttr) -> QueryResponse:
        """
        Perform a search query for AIA synoptic data.

        Parameters:
            query (SimpleAttr): The query parameters including time, wavelength, and cadence.

        Returns:
            QueryResponse: A response object containing the result of the query.
        """
        time_range = None
        wavelengths = []
        cadence_seconds = None

        for q in query:
            if isinstance(q, attrs.Time):
                time_range = q
            elif isinstance(q, attrs.Wavelength):
                wavelength_value = q.min.to(u.angstrom).value
                wavelengths.append(int(wavelength_value))
            elif isinstance(q, attrs.Sample):
                cadence_seconds = q.value

        if not time_range:
            logger.error("Time range must be specified for the AIASynopticClient.")
            raise ValueError("Time range must be specified for the AIASynopticClient.")

        if not wavelengths:
            wavelengths = self.known_wavelengths
        else:
            wavelengths = [str(wl).zfill(4) for wl in wavelengths]

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
            "Instrument": [],
            "Wavelength": [],
            "url": [],
            "fileid": [],
        }
        for url in urls:
            filename = os.path.basename(url)
            name_part = filename[3:-5]
            parts = name_part.split("_")
            if len(parts) == 3:
                date_str = parts[0]
                time_str = parts[1]
                wavelength_str = parts[2]
                datetime_str = date_str + time_str
                try:
                    start_time = datetime.strptime(datetime_str, "%Y%m%d%H%M")
                except ValueError:
                    start_time = None
                data["Start Time"].append(start_time)
                data["Wavelength"].append(int(wavelength_str))
            else:
                data["Start Time"].append(None)
                data["Wavelength"].append(None)
            data["Instrument"].append("AIA")
            data["url"].append(url)
            data["fileid"].append(filename)
        table = QueryResponseTable(data, client=self)
        return QueryResponse(table)

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
        current_time = time_range.start.datetime
        end_time = time_range.end.datetime
        urls = []

        if cadence_seconds is not None:
            cadence_timedelta = timedelta(seconds=cadence_seconds)
        else:
            cadence_timedelta = timedelta(minutes=1)

        while current_time <= end_time:
            for wavelength in wavelengths:
                formatted_baseurl = current_time.strftime(self.baseurl)
                formatted_url = current_time.strftime(self.pattern).format(
                    baseurl=formatted_baseurl, wavelength=str(wavelength).zfill(4)
                )
                urls.append(formatted_url)
            current_time += cadence_timedelta

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
            max_conn = kwargs.get("max_conn", 10)
            downloader.max_conn = max_conn
            for record in query_result:
                downloader.enqueue_file(record["url"], path=download_path)
            return downloader.download()
        except Exception as e:
            logger.error(f"Error while fetching data: {e}")
            raise


# Register the client with Fido
if AIASynopticClient not in Fido.registry:
    Fido.registry[AIASynopticClient] = AIASynopticClient._can_handle_query
    logger.info("Synoptic Fido Client Loaded!")
