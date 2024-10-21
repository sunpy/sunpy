from sunpy.net.dataretriever import GenericClient
from sunpy.net import attrs as a
from sunpy.net.scraper import Scraper
import astropy.units as u
from datetime import timedelta
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

    baseurl = "https://jsoc1.stanford.edu/data/aia/synoptic/%Y/%m/%d/H%H00/AIA%Y%m%d_%H%M_{wavelength:04d}.fits"
    pattern = r".*/aia/synoptic/\d{4}/\d{2}/\d{2}/H\d{4}/AIA\d{8}_\d{4}_(\d{4}).fits"
    known_wavelengths = [94, 131, 171, 193, 211, 304, 335, 1600, 1700, 4500]

    required = {a.Time, a.Instrument}
    optional = {a.Sample, a.Level, a.Wavelength, a.ExtentType}

    @property
    def info_url(self):
        return "https://jsoc1.stanford.edu/data/aia/synoptic/"

    @classmethod
    def register_values(cls):
        adict = {
            a.Instrument: [("AIA", "Data from the Atmospheric Imaging Assembly instrument.")],
            a.Physobs: [("intensity", "Brightness or intensity of the solar atmosphere at different wavelengths.")],
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
            print("BROKE1")
            return False

        # Ensure there are no unsupported attributes in the query
        if not query_attrs.issubset(supported_attrs):
            print("BROKE2")
            return False

        return True

    def search(self, *args, **kwargs):
        """
        Perform a search for AIA synoptic data.
        """
        matchdict = self._get_match_dict(*args, **kwargs)
        time = matchdict["time"]
        wavelength = matchdict.get("wavelength", None)

        # If no wavelength is provided, search for all known wavelengths
        wavelengths_to_search = [wavelength] if wavelength else self.known_wavelengths

        all_results = []
        for wl in wavelengths_to_search:
            self.scraper = self._get_scraper(wl)
            result = self.scraper.filelist(time)
            all_results.extend(result)

        # If there are no results, return an empty QueryResponse
        if not all_results:
            return QueryResponse([], client=self)

        # Convert file list to QueryResponse
        query_response = self._make_records(all_results, matchdict)

        # Remove the 'URL' column from the QueryResponse
        query_response.remove_column("URL")

        return query_response

    def _get_scraper(self, wavelength):
        """
        Return a Scraper instance for the given wavelength.
        """
        url = self.baseurl.format(wavelength=wavelength)
        return Scraper(url, regex=self.pattern)

    def _make_records(self, filelist, matchdict):
        """
        Convert a list of file URLs to a QueryResponse object.
        """
        records = []
        last_time = None  # Track the time of the last appended record
        sample = matchdict.get("sample", None)  # Get the sampling interval if provided

        # Convert sample interval to minutes if provided
        sample_minutes = sample.value/60 if sample else None

        for file in filelist:
            time_match = self.scraper._extract_date(file)
            if time_match:
                start_time = time_match.datetime
                # Apply sampling logic: skip if the sample interval is not met
                if sample_minutes is not None:
                    if last_time is not None and (start_time - last_time).total_seconds() / 60 < sample_minutes:
                        continue  # Skip appending this record

                # Update the last_time to the current start_time if appending
                last_time = start_time

                # Extract metadata from matchdict if available
                wavelength = matchdict.get("wavelength", None)
                # Create a record for each file
                record = {
                    'Start Time': start_time,
                    'End Time': start_time + timedelta(seconds=119),  # Add 1 minute and 59 seconds
                    'Instrument': 'AIA',
                    'Physobs': 'intensity',
                    'Source': 'SDO',
                    'Provider': 'JSOC',
                    'Level': 'synoptic',
                    'ExtentType': 'synoptic',
                    'Wavelength': f"{wavelength} Å" if wavelength else 'Unknown',
                    'URL': file
                }
                records.append(record)

        # Convert the list of records to a QueryResponse
        return QueryResponse(records, client=self)

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
                matchdict["sample"] = arg

        return matchdict