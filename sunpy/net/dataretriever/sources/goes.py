# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014

import os
from datetime import datetime
from itertools import compress
from urllib.parse import urlsplit

import astropy.units as u
from astropy.time import Time, TimeDelta

from sunpy import config
from sunpy.net.dataretriever import GenericClient
from sunpy.time import TimeRange, parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring
from sunpy.util.scraper import Scraper

TIME_FORMAT = config.get("general", "time_format")

__all__ = ["XRSClient", "SUVIClient"]


class XRSClient(GenericClient):
    """
    Provides access to the GOES XRS fits files archive.

    Searches data hosted by the `Solar Data Analysis Center <https://umbra.nascom.nasa.gov/goes/fits/>`__.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.xrs)  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the XRSClient:
         Start Time           End Time      Source Instrument Wavelength
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
        timerange : `~sunpy.time.TimeRange`
            The time range you want the files for.

        Returns
        -------
        `list`
            The URL(s) for the corresponding timerange.
        """
        timerange = TimeRange(timerange.start.strftime('%Y-%m-%d'), timerange.end)
        if timerange.end < parse_time("1999/01/15"):
            goes_file = "%Y/go{satellitenumber:02d}%y%m%d.fits"
        elif timerange.start < parse_time("1999/01/15") and timerange.end >= parse_time("1999/01/15"):
            return self._get_overlap_urls(timerange)
        else:
            goes_file = "%Y/go{satellitenumber}%Y%m%d.fits"

        goes_pattern = f"https://umbra.nascom.nasa.gov/goes/fits/{goes_file}"
        satellitenumber = kwargs.get("satellitenumber", self._get_goes_sat_num(timerange.start))
        goes_files = Scraper(goes_pattern, satellitenumber=satellitenumber)

        return goes_files.filelist(timerange)

    def _get_overlap_urls(self, timerange):
        """
        Return a list of URLs over timerange when the URL path changed format `%Y` to `%y`
        on the date 1999/01/15

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            The time range you want the files for.
        Returns
        -------
        `list`
            The URL(s) for the corresponding timerange.
        """
        tr_before = TimeRange(timerange.start, parse_time("1999/01/14"))
        tr_after = TimeRange(parse_time("1999/01/15"), timerange.end)
        urls_before = self._get_url_for_timerange(tr_before)
        urls_after = self._get_url_for_timerange(tr_after)
        return urls_before + urls_after

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

    @classmethod
    def _attrs_module(cls):
        return 'goes', 'sunpy.net.dataretriever.attrs.goes'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        goes_number = [2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        adict = {attrs.Instrument: [
            ("GOES", "The Geostationary Operational Environmental Satellite Program."),
            ("XRS", "GOES X-ray Flux")],
            attrs.goes.SatelliteNumber: [(str(x), f"GOES Satellite Number {x}") for x in goes_number]}
        return adict


class SUVIClient(GenericClient):
    """
    Provides access to data from the GOES Solar Ultraviolet Imager (SUVI).

    SUVI data are provided by NOAA at the following url
    https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/
    The SUVI instrument was first included on GOES-16. It produces level-1b as
    well as level-2 data products. Level-2 data products are a weighted average
    of level-1b product files and therefore provide higher imaging dynamic
    range than individual images. The exposure time of level 1b images range
    from 1 s to 0.005 s. SUVI supports the following wavelengths;
    94, 131, 171, 195, 284, 304 angstrom. If no wavelength is specified, images
    from all wavelengths are returned.

    Note
    ----
    GOES-16 began providing regular level-1b data on 2018-06-01.  At the time
    of writing, SUVI on GOES-17 is operational but currently does not provide
    Level-2 data.
    """

    @add_common_docstring(**_variables_for_parse_time_docstring())
    def _get_goes_sat_num(self, date):
        """
        Determines the best satellite number for a given date.

        Parameters
        ----------
        date : {parse_time_types}
            The date to determine which satellite is active.

        Note
        ----
        At the time this function was written.
        GOES-17 is operational but currently does not provide Level 2 data therefore it is never returned.
        The GOES-16 start date is based on the availability of regular level 1b data.
        """

        # GOES-17 is operational but currently does not provide Level 2 data
        # GOES-16 start date is based on the availability of regular level 1b data
        suvi_operational = {
            16: TimeRange("2018-06-01", parse_time("now")),
        }

        results = []
        for sat_num in suvi_operational:
            if date in suvi_operational[sat_num]:
                # if true then the satellite with sat_num is available
                results.append(sat_num)

        if results:
            # Return the newest satellite
            return max(results)
        else:
            # if no satellites were found then raise an exception
            raise ValueError(f"No operational SUVI instrument on {date.strftime(TIME_FORMAT)}")

    def _get_time_for_url(self, urls):
        these_timeranges = []

        for this_url in urls:
            if this_url.count('/l2/') > 0:  # this is a level 2 data file
                start_time = parse_time(os.path.basename(this_url).split('_s')[2].split('Z')[0])
                end_time = parse_time(os.path.basename(this_url).split('_e')[1].split('Z')[0])
                these_timeranges.append(TimeRange(start_time, end_time))
            if this_url.count('/l1b/') > 0:  # this is a level 1b data file
                start_time = datetime.strptime(os.path.basename(this_url).split('_s')[
                                               1].split('_e')[0][:-1], '%Y%j%H%M%S')
                end_time = datetime.strptime(os.path.basename(this_url).split('_e')[
                                             1].split('_c')[0][:-1], '%Y%j%H%M%S')
                these_timeranges.append(TimeRange(start_time, end_time))
        return these_timeranges

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns urls to the SUVI data for the specified time range.

        Parameters
        ----------
        timerange: `sunpy.time.TimeRange`
            Time range for which data is to be downloaded.
        level : `str`, optional
            The level of the data. Possible values are 1b and 2 (default).
        wavelength : `astropy.units.Quantity` or `tuple`, optional
            Wavelength band. If not given, all wavelengths are returned.
        satellitenumber : `int`, optional
            GOES satellite number. Must be >= 16. Default is 16.
        """
        base_url = "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes{goes_number}/"
        supported_waves = [94, 131, 171, 195, 284, 304]
        supported_levels = ("2", "1b")

        # these are optional requirements so if not provided assume defaults
        # if wavelength is not provided assuming all of them
        if "wavelength" in kwargs.keys():
            wavelength_input = kwargs.get("wavelength")
            if isinstance(wavelength_input, u.Quantity):  # not a range
                if int(wavelength_input.to_value('Angstrom')) not in supported_waves:
                    raise ValueError(f"Wavelength {kwargs.get('wavelength')} not supported.")
                else:
                    wavelength = [kwargs.get("wavelength")]
            else:  # Range was provided
                compress_index = [wavelength_input.wavemin <= this_wave <=
                                  wavelength_input.wavemax for this_wave in (supported_waves * u.Angstrom)]
                if not any(compress_index):
                    raise ValueError(
                        f"Wavelength {wavelength_input} not supported.")
                else:
                    wavelength = list(compress(supported_waves, compress_index)) * u.Angstrom
        else:  # no wavelength provided return all of them
            wavelength = supported_waves * u.Angstrom
        # check that the input wavelength can be converted to angstrom
        waves = [int(this_wave.to_value('angstrom', equivalencies=u.spectral()))
                 for this_wave in wavelength]
        # use the given satellite number or choose the best one
        satellitenumber = int(kwargs.get(
            "satellitenumber", self._get_goes_sat_num(timerange.start)))
        if satellitenumber < 16:
            raise ValueError(f"Satellite number {satellitenumber} not supported.")
        # default to the highest level of data
        level = str(kwargs.get("level", "2"))  # make string in case the input is a number

        if level not in supported_levels:
            raise ValueError(f"Level {level} is not supported.")

        results = []
        for this_wave in waves:
            if level == "2":
                search_pattern = base_url + \
                    r'l{level}/data/suvi-l{level}-ci{wave:03}/%Y/%m/%d/dr_suvi-l{level}-ci{wave:03}_g{goes_number}_s%Y%m%dT%H%M%SZ_.*\.fits'
            elif level == "1b":
                if this_wave in [131, 171, 195, 284]:
                    search_pattern = base_url + \
                        r'l{level}/suvi-l{level}-fe{wave:03}/%Y/%m/%d/OR_SUVI-L{level}-Fe{wave:03}_G{goes_number}_s%Y%j%H%M%S.*\.fits.gz'
                elif this_wave == 304:
                    search_pattern = base_url + \
                        r'l{level}/suvi-l{level}-he{wave:03}/%Y/%m/%d/OR_SUVI-L{level}-He{wave_minus1:03}_G{goes_number}_s%Y%j%H%M%S.*\.fits.gz'
                elif this_wave == 94:
                    search_pattern = base_url + \
                        r'l{level}/suvi-l{level}-fe{wave:03}/%Y/%m/%d/OR_SUVI-L{level}-Fe{wave_minus1:03}_G{goes_number}_s%Y%j%H%M%S.*\.fits.gz'

            if search_pattern.count('wave_minus1'):
                scraper = Scraper(search_pattern, level=level, wave=this_wave,
                                  goes_number=satellitenumber, wave_minus1=this_wave-1)
            else:
                scraper = Scraper(search_pattern, level=level, wave=this_wave,
                                  goes_number=satellitenumber)
            results.extend(scraper.filelist(timerange))
        return results

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
            All specified query objects.

        Returns
        -------
        `bool`
            answer as to whether client can service the query.
        """
        # Import here to prevent circular imports
        from sunpy.net import attrs as a

        required = {a.Time, a.Instrument}
        optional = {a.Wavelength, a.Level, a.goes.SatelliteNumber}
        all_attrs = {type(x) for x in query}

        ops = all_attrs - required
        # check to ensure that all optional requirements are in approved list
        if ops and not all(elem in optional for elem in ops):
            return False

        # if we get this far we have either Instrument and Time
        # or Instrument, Time and Wavelength
        check_var_count = 0
        for x in query:
            if isinstance(x, a.Instrument) and x.value.lower() == 'suvi':
                check_var_count += 1

        if check_var_count == 1:
            return True
        else:
            return False

    @classmethod
    def _attrs_module(cls):
        return 'goes', 'sunpy.net.dataretriever.attrs.goes'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        goes_number = [16, 17]
        adict = {attrs.Instrument: [
            ("SUVI", "The Geostationary Operational Environmental Satellite Program.")],
            attrs.goes.SatelliteNumber: [(str(x), f"GOES Satellite Number {x}") for x in goes_number]}
        return adict
