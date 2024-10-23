import numpy as np
import re
from datetime import datetime, timedelta
from sunpy.net.dataretriever import GenericClient
from sunpy.net import attrs as a
import astropy.units as u
from sunpy.net.dataretriever.client import QueryResponse

__all__ = ["AIASynopticClient"]


class AIASynopticClient(GenericClient):
    """
    A client for retrieving AIA synoptic data from JSOC.

    This client retrieves synoptic AIA data from the following source:
    https://jsoc1.stanford.edu/data/aia/synoptic/

    The synoptic dataset includes lower-resolution (1k) images, with additional
    image processing steps like downsampling and time integration.
    """

    baseurl = "https://jsoc1.stanford.edu/data/aia/synoptic/"
    known_wavelengths = [94, 131, 171, 193, 211, 304, 335, 1600, 1700, 4500]

    required = {a.Time, a.Instrument, a.Level}
    supported = required + {a.Sample, a.Wavelength, a.ExtentType}

    @property
    def info_url(self):
        return self.baseurl

    @classmethod
    def register_values(cls):
        adict = {
            a.Instrument: [("AIA", "Data from the Atmospheric Imaging Assembly instrument.")],
            a.Physobs: [
                ("intensity", "Brightness or intensity of the solar atmosphere at different wavelengths.")
            ],
            a.Source: [("SDO", "The Solar Dynamics Observatory.")],
            a.Wavelength: [(f"{wv:04d}", f"{wv} Å") for wv in cls.known_wavelengths],
            a.Provider: [("JSOC", "Joint Science Operations Center at Stanford.")],
            a.Level: [("synoptic", "Level 1.5 data processed for quicker analysis.")],
            a.ExtentType: [("synoptic", "Level 1.5 data processed for quicker analysis.")],
        }
        return adict

    def _can_handle_query(self, *query):
        """
        Determines whether this client can handle the given query.
        """
        required_attrs = {a.Instrument, a.Level}
        supported_attrs = {a.Time, a.Instrument, a.Sample, a.Level, a.Wavelength, a.ExtentType}

        # Check if all required attributes are present in the query
        query_attrs = set(type(attr) for attr in query)
        if not required_attrs.issubset(query_attrs):
            return False

        # Ensure there are no unsupported attributes in the query
        if not query_attrs.issubset(supported_attrs):
            return False

        return True

    def search(self, *args, **kwargs):
        """
        Perform a search for AIA synoptic data.
        """
        matchdict = self._get_match_dict(*args, **kwargs)
        time_range = matchdict["time"]
        start_time = time_range.start.datetime
        end_time = time_range.end.datetime
        sample = matchdict.get("sample", None)  # Get the sampling interval if provided
        sample_minutes = (
            sample.to(u.min).value if sample is not None else 2
        )  # Default to 2-minute intervals if not provided

        # List of wavelengths to search
        wavelengths = [matchdict["wavelength"]] if "wavelength" in matchdict else self.known_wavelengths

        # Calculate the number of steps based on the sample interval
        total_minutes = (end_time - start_time).total_seconds() / 60
        num_steps = int(np.ceil(total_minutes / sample_minutes))

        # Generate the array of times
        time_steps = [
            start_time + timedelta(minutes=i * sample_minutes)
            for i in range(num_steps + 1)
            if start_time + timedelta(minutes=i * sample_minutes) <= end_time
        ]

        # Generate all possible file URLs based on the time steps and wavelengths
        all_results = []
        for current_time in time_steps:
            for wl in wavelengths:
                # Format the URL based on the known structure
                url = (
                    f"{self.baseurl}"
                    f"{current_time.year}/{current_time:%m}/{current_time:%d}/H{current_time:%H}00/"
                    f"AIA{current_time:%Y%m%d}_{current_time:%H%M}_{wl:04d}.fits"
                )
                all_results.append(url)

        # If there are no results, return an empty QueryResponse
        if not all_results:
            return QueryResponse([], client=self)

        # Convert file list to QueryResponse
        query_response = self._make_records(all_results, matchdict)

        # Remove the 'URL' column from the QueryResponse
        # query_response.remove_column("URL")

        return query_response

    def _make_records(self, all_results, matchdict):
        """
        Convert a list of file URLs to a QueryResponse object.
        """
        records = []
        for url in all_results:
            # Extract the timestamp from the URL
            time_match = self._extract_date_from_url(url)

            if time_match:
                start_time = time_match
                wavelength = self._extract_wavelength(url)
                # Extract metadata from matchdict if available
                # wavelength = matchdict.get("wavelength", None)
                # Create a record for each file
                record = {
                    "Start Time": start_time,
                    "End Time": start_time + timedelta(seconds=119),  # Add 1 minute and 59 seconds
                    "Instrument": "AIA",
                    "Physobs": "intensity",
                    "Source": "SDO",
                    "Provider": "JSOC",
                    "Level": "synoptic",
                    "ExtentType": "synoptic",
                    "Wavelength": f"{wavelength} Å" if wavelength else "Unknown",
                    "url": url,  # Make sure 'url' is lowercase
                }
                records.append(record)

        # Convert the list of records to a QueryResponse
        return QueryResponse(records, client=self)

    def _extract_date_from_url(self, url):
        """
        Extract the date and time from the given URL based on the expected pattern.
        """
        pattern = r".*/(\d{4})/(\d{2})/(\d{2})/H(\d{2})00/AIA(\d{8})_(\d{4})_\d{4}.fits"
        match = re.search(pattern, url)
        if match:
            year, month, day, hour, yyyymmdd, hhmm = match.groups()
            return datetime.strptime(f"{yyyymmdd}_{hhmm}", "%Y%m%d_%H%M")
        return None

    def _extract_wavelength(self, url: str) -> int | None:
        """
        Extract the 4-digit wavelength value from the given URL based on the expected pattern,
        set it to self.wavelength, and return the wavelength as an integer.

        Args:
            url (str): The URL string from which to extract the wavelength.

        Returns:
            int | None: The extracted wavelength as an integer if successful, otherwise None.
        """
        next = url.split("_")[-1].split(".")[0]
        return next if next else None

    def _get_match_dict(self, *args, **kwargs):
        """
        Extract parameters from the query for use in setting up the scraper.
        """
        matchdict = {}
        for arg in args:
            if isinstance(arg, a.Time):
                matchdict["time"] = arg
            if isinstance(arg, a.Wavelength):
                wavelength_value = int(arg.min.to(u.angstrom).value)
                if wavelength_value in self.known_wavelengths:
                    matchdict["wavelength"] = wavelength_value
                else:
                    raise ValueError(f"Wavelength {wavelength_value} Å is not a known wavelength for AIA.")
            if isinstance(arg, a.Sample):
                # Ensure that the sample is stored as a Quantity with time units
                if isinstance(arg.value, u.Quantity):
                    matchdict["sample"] = arg.value
                else:
                    matchdict["sample"] = arg.value * u.s  # Default to seconds if units are missing

        return matchdict
