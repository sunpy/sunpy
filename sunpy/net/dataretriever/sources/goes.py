# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014

import os
from datetime import datetime
from itertools import compress

import astropy.units as u

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
    4 Results from the XRSClient:
         Start Time           End Time      Instrument ... Provider SatelliteNumber
    ------------------- ------------------- ---------- ... -------- ---------------
    2016-01-01 00:00:00 2016-01-01 23:59:59        XRS ...     SDAC              13
    2016-01-02 00:00:00 2016-01-02 23:59:59        XRS ...     SDAC              13
    2016-01-01 00:00:00 2016-01-01 23:59:59        XRS ...     SDAC              15
    2016-01-02 00:00:00 2016-01-02 23:59:59        XRS ...     SDAC              15
    <BLANKLINE>
    <BLANKLINE>

    """
    baseurl = r'https://umbra.nascom.nasa.gov/goes/fits/%Y/go(\d){2}(\d){6,8}\.fits'
    pattern = '{}/fits/{year:4d}/go{SatelliteNumber:02d}{}{month:2d}{day:2d}.fits'

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
            attrs.goes.SatelliteNumber: [(str(x), f"GOES Satellite Number {x}") for x in goes_number],
            attrs.Source: [('NASA', 'The National Aeronautics and Space Administration.')],
            attrs.Physobs: [('irradiance', 'the flux of radiant energy per unit area.')],
            attrs.Provider: [('SDAC', 'The Solar Data Analysis Center.')]}
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

    @classmethod
    def _attrs_module(cls):
        return 'goes', 'sunpy.net.dataretriever.attrs.goes'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        goes_number = [16, 17]
        adict = {attrs.Instrument: [
            ("SUVI", "GOES Solar Ultraviolet Imager.")],
            attrs.goes.SatelliteNumber: [(str(x), f"GOES Satellite Number {x}") for x in goes_number],
            attrs.Source: [('GOES', 'The Geostationary Operational Environmental Satellite Program.')],
            attrs.Physobs: [('flux', 'a measure of the amount of radiation received by an object from a given source.')],
            attrs.Provider: [('NOAA', 'The National Oceanic and Atmospheric Administration.')]}
        return adict
