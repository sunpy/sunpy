# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014

import os
from urllib.parse import urlsplit
from datetime import timedelta
import warnings

from astropy.time import TimeDelta
from astropy.time import Time
import astropy.units as u

from sunpy.time import parse_time, TimeRange
from sunpy import config
from sunpy.util.scraper import Scraper
from ..client import GenericClient
from sunpy import log
from sunpy.util.exceptions import SunpyUserWarning
from sunpy.net import attrs as a


TIME_FORMAT = config.get("general", "time_format")

__all__ = ["XRSClient", "SUVIClient"]


class XRSClient(GenericClient):
    """
    Provides access to the GOES XRS fits files
    `archive <https://umbra.nascom.nasa.gov/goes/fits/>`__ hosted
    by the `Solar Data Analysis Center <https://umbra.nascom.nasa.gov/index.html/>`__.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument('XRS'))  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the XRSClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str4     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2016-01-01 00:00:00 2016-01-01 23:59:59   nasa       goes        nan
    2016-01-02 00:00:00 2016-01-02 23:59:59   nasa       goes        nan
    <BLANKLINE>
    <BLANKLINE>

    """

    def _get_goes_sat_num(self, date):
        """
        Determines the satellite number for a given date.

        Parameters
        ----------

        date : `astropy.time.Time`
            The date to determine which satellite is active.
        """
        goes_operational = {
            2: TimeRange("1981-01-01", "1983-04-30"),
            5: TimeRange("1983-05-02", "1984-07-31"),
            6: TimeRange("1983-06-01", "1994-08-18"),
            7: TimeRange("1994-01-01", "1996-08-13"),
            8: TimeRange("1996-03-21", "2003-06-18"),
            9: TimeRange("1997-01-01", "1998-09-08"),
            10: TimeRange("1998-07-10", "2009-12-01"),
            11: TimeRange("2006-06-20", "2008-02-15"),
            12: TimeRange("2002-12-13", "2007-05-08"),
            13: TimeRange("2006-08-01", "2006-08-01"),
            14: TimeRange("2009-12-02", "2010-10-04"),
            15: TimeRange("2010-09-01", parse_time("now")),
        }

        results = []
        for sat_num in goes_operational:
            if date in goes_operational[sat_num]:
                # if true then the satellite with sat_num is available
                results.append(sat_num)

        if results:
            # Return the newest satellite
            return max(results)
        else:
            # if no satellites were found then raise an exception
            raise ValueError(
                "No operational GOES satellites on {}".format(
                    date.strftime(TIME_FORMAT)
                )
            )

    def _get_time_for_url(self, urls):
        times = []
        for uri in urls:
            uripath = urlsplit(uri).path

            # Extract the yymmdd or yyyymmdd timestamp
            datestamp = os.path.splitext(os.path.split(uripath)[1])[0][4:]

            # 1999-01-15 as an integer.
            if int(datestamp) <= 990115:
                start = Time.strptime(datestamp, "%y%m%d")
            else:
                start = Time.strptime(datestamp, "%Y%m%d")

            almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
            times.append(TimeRange(start, start + almost_day))

        return times

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns a URL to the GOES data for the specified date.

        Parameters
        ----------
        timerange: sunpy.time.TimeRange
            time range for which data is to be downloaded.
        satellitenumber : int
            GOES satellite number (default = 15)
        data_type : str
            Data type to return for the particular GOES satellite. Supported
            types depend on the satellite number specified. (default = xrs_2s)
        """
        # find out which satellite and datatype to query from the query times
        base_url = "https://umbra.nascom.nasa.gov/goes/fits/"
        start_time = Time(timerange.start.strftime("%Y-%m-%d"))
        # make sure we are counting a day even if only a part of it is in the query range.
        day_range = TimeRange(
            timerange.start.strftime("%Y-%m-%d"), timerange.end.strftime("%Y-%m-%d")
        )
        total_days = int(day_range.days.value) + 1
        result = list()

        # Iterate over each day in the input timerange and generate a URL for
        # it.
        for day in range(total_days):
            # It is okay to convert to datetime here as the start_time is a date
            # hence we don't necesserily gain anything.
            # This is necessary because when adding a day to a Time, we may
            # end up with the same day if the day is a leap second day
            date = start_time.datetime + timedelta(days=day)
            regex = date.strftime("%Y") + "/go{sat:02d}"
            if date < parse_time("1999/01/15"):
                regex += date.strftime("%y%m%d") + ".fits"
            else:
                regex += date.strftime("%Y%m%d") + ".fits"
            satellitenumber = kwargs.get(
                "satellitenumber", self._get_goes_sat_num(date)
            )
            url = base_url + regex.format(sat=satellitenumber)
            result.append(url)
        return result

    def _makeimap(self):
        """
        Helper function used to hold information about source.
        """
        self.map_["source"] = "nasa"
        self.map_["instrument"] = "goes"
        self.map_["physobs"] = "irradiance"
        self.map_["provider"] = "sdac"

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Answers whether client can service the query.

        Parameters
        ----------
        query : list of query objects

        Returns
        -------
        boolean
            answer as to whether client can service the query
        """
        chkattr = ["Time", "Instrument", "SatelliteNumber"]
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == "Instrument" and x.value.lower() in (
                "xrs",
                "goes",
            ):
                return all(chklist)
        return False


class SUVIClient(GenericClient):
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns a URL to the SUVI data for the specified time range.

        Parameters
        ----------
        timerange: sunpy.time.TimeRange
            time range for which data is to be downloaded.
        level : str
            The level of the data.
        wavelength :
            Wavelength band.
        satellitenumber : int
            GOES satellite number. Must be >= 16.

        """
        base_url = "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes{goes_number}/"

        if not kwargs.get("wavelength", None):
            # should never reach here if called by dataretriever
            raise ValueError("Queries must specify a wavelength.")
        else:
            wavelength = kwargs["wavelength"]

        wavelength = wavelength.to(u.Angstrom, equivalencies=u.spectral())
        wave = int(wavelength.to_value('angstrom'))

        # these are optional requirements so if not provided assume defaults
        satellitenumber = int(kwargs.get("satellitenumber", 16))
        level = str(kwargs.get("level", "2"))  # cast to string because it may be provided as number

        if level == "2":
            search_pattern = base_url + 'l{level}/data/suvi-l{level}-ci{wave:03}/%Y/%m/%d/dr_suvi-l{level}-ci{wave:03}_g16_s%Y%m%dT%H%M%SZ_.*\.fits'
        elif level == "1b":
            if wave in [131, 171, 195, 284]:
                search_pattern = base_url + 'l{level}/suvi-l{level}-fe{wave:03}/%Y/%m/%d/OR_SUVI-L{level}-Fe{wave:03}_G16_s%Y%j%H%M%S.*\.fits.gz'
            elif wave == 304:
                search_pattern = base_url + 'l{level}/suvi-l{level}-he{wave:03}/%Y/%m/%d/OR_SUVI-L{level}-He{wave_minus1:03}_G16_s%Y%j%H%M%S.*\.fits.gz'
            elif wave == 94:
                search_pattern = base_url + 'l{level}/suvi-l{level}-fe{wave:03}/%Y/%m/%d/OR_SUVI-L{level}-Fe{wave_minus1:03}_G16_s%Y%j%H%M%S.*\.fits.gz'
        else:  # should never reach here if called by dataretriever
            return []

        if search_pattern.count('wave_minus1'):
            scraper = Scraper(search_pattern, level=level, wave=wave,
                              goes_number=satellitenumber, wave_minus1=wave-1)
        else:
            scraper = Scraper(search_pattern, level=level, wave=wave,
                              goes_number=satellitenumber)
        return scraper.filelist(timerange)

    def _makeimap(self):
        """
        Helper Function used to hold information about source.
        """
        self.map_['source'] = 'GOES'
        self.map_['provider'] = 'NOAA'
        self.map_['instrument'] = 'SUVI'
        self.map_['physobs'] = 'flux'

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Answers whether client can service the query.

        Parameters
        ----------
        query : `tuple`
            All specified query objects

        Returns
        -------
        `bool`
            answer as to whether client can service the query
        """
        supported_levels = ("2", "1b")
        supported_waves = [94, 131, 171, 195, 284, 304]
        minimum_supported_satellitenumber = 16

        required = {a.Time, a.Instrument, a.Wavelength}
        optional = {} #, a.SatelliteNumber, a.Level}
        all_attrs = {type(x) for x in query}

        ops = all_attrs - required
        # If ops is empty or equal to optional we are ok, otherwise we don't
        # match
        if ops and ops != optional:
            return False

        # if we get this far we have either Instrument and Time
        # or Instrument, Time and Wavelength
        check_var_count = 0
        for x in query:
            if isinstance(x, a.Instrument) and x.value.lower() == 'suvi':
                check_var_count += 1
            if isinstance(a, a.Wavelength) and int(x.value.to_value('Angstrom')) is in supported_waves:
                check_var_count += 1

            # how to check if level is appropriate, what is the appropriate attr?
            # how to check if the satellite number is appropriate, what is the appropriate attr?

        if check_var_count == 2:
            return True
        else
            return False
