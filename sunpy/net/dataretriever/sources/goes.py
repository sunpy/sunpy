# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014

from datetime import datetime
from collections import OrderedDict

import astropy.units as u
from astropy.time import Time

from sunpy import config
from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.net.scraper import Scraper
from sunpy.net.scraper_utils import get_timerange_from_exdict
from sunpy.time import TimeRange, parse_time

TIME_FORMAT = config.get("general", "time_format")

__all__ = ["XRSClient", "SUVIClient"]


class XRSClient(GenericClient):
    """
    Provides access to several GOES XRS files archive.

    For older GOES satellites (7 and below), NASA servers are used.
    For newer GOES satellites (8 and above), NOAA servers are used.

    For GOES 8-15, the data is the re-processed science-quality data.

    .. note::

        The new science quality data have the scaling factors removed for GOES 8-15
        and they are not added to GOES 16 AND 17 data products.
        This means the peak flux will be different to the older version of the data,
        such as those collected from the NASA servers.

    See the following readmes about these data:

    * GOES 1 - 7: https://umbra.nascom.nasa.gov/goes/fits/goes_fits_files_notes.txt

    * Reprocessed 8 - 15: https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/GOES_1-15_XRS_Science-Quality_Data_Readme.pdf

    * GOES-R 16 - 17: https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/docs/GOES-R_XRS_L2_Data_Readme.pdf

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.xrs)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    8 Results from the XRSClient:
    Source: <8: https://umbra.nascom.nasa.gov/goes/fits
    8-15: https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/
    16-17: https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/
    <BLANKLINE>
           Start Time               End Time        Instrument  Physobs   Source Provider Resolution SatelliteNumber filename_res
    ----------------------- ----------------------- ---------- ---------- ------ -------- ---------- --------------- ------------
    2016-01-01 00:00:00.000 2016-01-01 23:59:59.999        XRS irradiance   GOES     NOAA      irrad              13         gxrs
    2016-01-02 00:00:00.000 2016-01-02 23:59:59.999        XRS irradiance   GOES     NOAA      irrad              13         gxrs
    2016-01-01 00:00:00.000 2016-01-01 23:59:59.999        XRS irradiance   GOES     NOAA      avg1m              13         xrsf
    2016-01-02 00:00:00.000 2016-01-02 23:59:59.999        XRS irradiance   GOES     NOAA      avg1m              13         xrsf
    2016-01-01 00:00:00.000 2016-01-01 23:59:59.999        XRS irradiance   GOES     NOAA      irrad              15         gxrs
    2016-01-02 00:00:00.000 2016-01-02 23:59:59.999        XRS irradiance   GOES     NOAA      irrad              15         gxrs
    2016-01-01 00:00:00.000 2016-01-01 23:59:59.999        XRS irradiance   GOES     NOAA      avg1m              15         xrsf
    2016-01-02 00:00:00.000 2016-01-02 23:59:59.999        XRS irradiance   GOES     NOAA      avg1m              15         xrsf
    <BLANKLINE>
    <BLANKLINE>
    """
    # GOES XRS data from NASA servers up to GOES 7.
    pattern_old = 'https://umbra.nascom.nasa.gov/goes/fits/{{year:4d}}/go{{SatelliteNumber:2d}}{{}}{{month:2d}}{{day:2d}}.fits'
    # The reprocessed 8-15 data should be taken from NOAA.
    pattern_new = ("https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/goes{SatelliteNumber:02d}/{filename_res}-l2-{Resolution}_science"
                   "/{{year:4d}}/{{month:2d}}/sci_{filename_res}-l2-{Resolution}_g{SatelliteNumber:2d}_d{{year:4d}}{{month:2d}}{{day:2d}}_{{}}.nc")
    # GOES-R Series 16-17 XRS data from NOAA.
    pattern_r = ("https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes{SatelliteNumber:2d}/l2/data/xrsf-l2-{Resolution}_science"
                 "/{{year:4d}}/{{month:2d}}/sci_xrsf-l2-{Resolution}_g{SatelliteNumber:2d}_d{{year:4d}}{{month:2d}}{{day:2d}}_{{}}.nc")

    @property
    def info_url(self):
        return ("<8: https://umbra.nascom.nasa.gov/goes/fits \n"
                "8-15: https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/ \n"
                "16-17: https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/")

    def post_search_hook(self, i, matchdict):
        tr = get_timerange_from_exdict(i)
        rowdict = OrderedDict()
        rowdict['Start Time'] = tr.start
        rowdict['Start Time'].format = 'iso'
        rowdict['End Time'] = tr.end
        rowdict['End Time'].format = 'iso'
        rowdict["Instrument"] = matchdict["Instrument"][0].upper()
        rowdict["Physobs"] = matchdict["Physobs"][0]
        rowdict["url"] = i["url"]
        rowdict["Source"] = matchdict["Source"][0]
        if i["url"].endswith(".fits"):  # for older FITS files
            rowdict["Provider"] = matchdict["Provider"][0]
        else:  # only Resolution attrs for the netcdf files
            rowdict["Provider"] = matchdict["Provider"][1]
            if "avg1m" in i["url"]:
                rowdict["Resolution"] = "avg1m"
            elif ("flx1s" in i["url"]) or ("irrad" in i["url"]):
                rowdict["Resolution"] = "flx1s"
            else:
                raise RuntimeError("Could not parse resolution from URL")
        return rowdict

    def search(self, *args, **kwargs):
        matchdict = self._get_match_dict(*args, **kwargs)
        # this is for the case when the timerange overlaps with the provider change.
        if matchdict["Start Time"] < "2001-03-01" and matchdict["End Time"] >= "2001-03-01":
            matchdict_before, matchdict_after = matchdict.copy(), matchdict.copy()
            matchdict_after["Start Time"] = parse_time('2001-03-01')
            matchdict_before["End Time"] = parse_time('2001-03-01')
            metalist_before = self._get_metalist(matchdict_before)
            metalist_after = self._get_metalist(matchdict_after)
            metalist = metalist_before + metalist_after
        else:
            metalist = self._get_metalist(matchdict)
        return QueryResponse(metalist, client=self)

    def _get_metalist_fn(self, matchdict, pattern, **kwargs):
        """
        Function to help get list of OrderedDicts.
        """
        metalist = []
        scraper = Scraper(format=pattern, **kwargs)
        tr = TimeRange(matchdict["Start Time"], matchdict["End Time"])
        filemeta = scraper._extract_files_meta(tr, matcher=matchdict)
        for i in filemeta:
            rowdict = self.post_search_hook(i, matchdict)
            rowdict.update(kwargs)
            metalist.append(rowdict)
        return metalist

    def _get_metalist(self, matchdict):
        """
        Function to get the list of OrderDicts.
        This makes it easier for when searching for overlapping providers.
        """
        metalist = []
        # The data before the re-processed GOES 8-15 data.
        if (matchdict["End Time"] < "2001-03-01") or (matchdict["End Time"] >= "2001-03-01" and matchdict["Provider"] == ["sdac"]):
            metalist += self._get_metalist_fn(matchdict, self.pattern_old)
        # New data from NOAA. It searches for both the high cadence and 1 minute average data.
        else:
            if matchdict["End Time"] >= "2017-02-07":
                for sat in matchdict["SatelliteNumber"]:
                    if int(sat) >= 16:  # here check for GOES 16 and 17
                        for res in matchdict["Resolution"]:
                            metalist += self._get_metalist_fn(matchdict, self.pattern_r, SatelliteNumber=int(sat), Resolution=res)

            if matchdict["End Time"] <= "2020-03-04":
                for sat in matchdict["SatelliteNumber"]:
                    if (int(sat) >= 8) & (int(sat) <= 15):  # here check for GOES 8-15
                        # The 1 minute average data is at a different base URL to that of the high cadence data which is why things are done this way.
                        for res in matchdict["Resolution"]:
                            if res == "avg1m":
                                filename_res = "xrsf"
                                resolution = "avg1m"
                                metalist += self._get_metalist_fn(matchdict, self.pattern_new, SatelliteNumber=int(sat),
                                                                  filename_res=filename_res, Resolution=resolution)
                            elif res == "flx1s":
                                filename_res = "gxrs"
                                resolution = "irrad"
                                metalist += self._get_metalist_fn(matchdict, self.pattern_new, SatelliteNumber=int(sat),
                                                                  filename_res=filename_res, Resolution=resolution)
                            else:
                                raise RuntimeError(f"{res}` is not an accepted resolution attrs for the XRSClient")
        return metalist

    @classmethod
    def _attrs_module(cls):
        return 'goes', 'sunpy.net.dataretriever.attrs.goes'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        goes_number = [2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        adict = {attrs.Instrument: [
            ("GOES", "The Geostationary Operational Environmental Satellite Program."),
            ("XRS", "GOES X-ray Sensor")],
            attrs.Physobs: [('irradiance', 'the flux of radiant energy per unit area.')],
            attrs.Source: [("GOES", "The Geostationary Operational Environmental Satellite Program.")],
            attrs.Provider: [('SDAC', 'The Solar Data Analysis Center.'),
                             ('NOAA', 'The National Oceanic and Atmospheric Administration.')],
            attrs.goes.SatelliteNumber: [(str(x), f"GOES Satellite Number {x}") for x in goes_number],
            attrs.Resolution: [('flx1s', 'High-cadence measurements XRS observation, 1s for GOES-R, 2s for GOES 13-15, 3s for GOES<13'),
                               ('avg1m', '1-minute averages of XRS measurements')]}
        return adict


class SUVIClient(GenericClient):
    """
    Provides access to data from the GOES Solar Ultraviolet Imager (SUVI).

    `SUVI data are provided by NOAA <https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/>`__.

    The SUVI instrument was first included on GOES-16. It produces level 1b as
    well as level-2 data products. Level 2 data products are a weighted average
    of level 1b product files and therefore provide higher imaging dynamic
    range than individual images. The exposure time of level 1b images range
    from 1 s to 0.005 s.

    SUVI supports the following wavelengths; 94, 131, 171, 195, 284, 304 angstrom.
    If no wavelength is specified, images from all wavelengths are returned.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> import astropy.units as u
    >>> results = Fido.search(a.Time("2020/7/10", "2020/7/10 00:10"), a.Instrument('suvi'),a.Level.two,
    ...                       a.goes.SatelliteNumber(16), a.Wavelength(304*u.Angstrom))  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the SUVIClient:
    Source: https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes
    <BLANKLINE>
           Start Time               End Time        Instrument Physobs Source Provider SatelliteNumber Level Wavelength
                                                                                                              Angstrom
    ----------------------- ----------------------- ---------- ------- ------ -------- --------------- ----- ----------
    2020-07-10 00:00:00.000 2020-07-10 00:04:00.000       SUVI    flux   GOES     NOAA              16     2      304.0
    2020-07-10 00:04:00.000 2020-07-10 00:08:00.000       SUVI    flux   GOES     NOAA              16     2      304.0
    2020-07-10 00:08:00.000 2020-07-10 00:12:00.000       SUVI    flux   GOES     NOAA              16     2      304.0
    <BLANKLINE>
    <BLANKLINE>
    """

    pattern1b = ('https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/'
                 'goes{SatelliteNumber}/l1b/suvi-l1b-{elem:02}{wave:03}/'
                 '{{year:4d}}/{{month:2d}}/{{day:2d}}/{{}}_s{{:7d}}{{hour:2d}}{{minute:2d}}{{second:2d}}'
                 '{{:1d}}_e{{:7d}}{{ehour:2d}}{{eminute:2d}}{{esecond:2d}}{{:1d}}_{{}}.fits.gz')
    pattern2 = ('https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/'
                'goes/goes{SatelliteNumber}/l2/data/suvi-l2-ci{wave:03}/{{year:4d}}/'
                '{{month:2d}}/{{day:2d}}/dr_suvi-l2-ci{{Wavelength:03d}}_g{SatelliteNumber}_s'
                '{{year:4d}}{{month:2d}}{{day:2d}}T{{hour:2d}}{{minute:2d}}{{second:2d}}Z_e'
                '{{eyear:4d}}{{emonth:2d}}{{eday:2d}}T{{ehour:2d}}{{eminute:2d}}{{esecond:2d}}Z_{{}}.fits')

    @property
    def info_url(self):
        return 'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes'

    def post_search_hook(self, i, matchdict):
        # Extracting start times and end times
        start = Time(datetime(i['year'], i['month'], i['day'], i['hour'], i['minute'], i['second']))
        start.format = 'iso'
        end = Time(datetime(i['year'], i['month'], i['day'], i['ehour'], i['eminute'], i['esecond']))
        end.format = 'iso'

        rowdict = OrderedDict()
        rowdict['Start Time'] = start
        rowdict['End Time'] = end
        rowdict['Instrument'] = matchdict['Instrument'][0].upper()
        rowdict['Physobs'] = matchdict['Physobs'][0]
        rowdict['Source'] = matchdict['Source'][0]
        rowdict['Provider'] = matchdict['Provider'][0]
        rowdict['url'] = i['url']
        return rowdict

    def search(self, *args, **kwargs):
        supported_waves = [94, 131, 171, 195, 284, 304]*u.Angstrom
        all_waves = []
        matchdict = self._get_match_dict(*args, **kwargs)
        req_wave = matchdict.get('Wavelength', None)
        if req_wave is not None:
            wmin_round = round(req_wave.min.value) * u.Angstrom
            wmax_round = round(req_wave.max.value) * u.Angstrom
            req_wave = a.Wavelength(wmin_round, wmax_round)
            for wave in supported_waves:
                if wave in req_wave:
                    all_waves.append(int(wave.value))
        else:
            all_waves = [int(i.value) for i in supported_waves]
        all_satnos = matchdict.get('SatelliteNumber')
        all_levels = matchdict.get('Level')
        metalist = []

        # iterating over all possible Attr values through loops
        for satno in all_satnos:
            for level in all_levels:
                for wave in all_waves:
                    formdict = {'wave': wave, 'SatelliteNumber': satno}
                    if str(level) == '1b':
                        formdict['elem'] = 'fe'
                        if wave == 304:
                            formdict['elem'] = 'he'
                        pattern = self.pattern1b
                    elif str(level) == '2':
                        pattern = self.pattern2
                    else:
                        raise ValueError(f"Level {level} is not supported.")
                    # formatting pattern using Level, SatelliteNumber and Wavelength
                    urlpattern = pattern.format(**formdict)
                    urlpattern = urlpattern.replace('{', '{{').replace('}', '}}')
                    scraper = Scraper(format=urlpattern)
                    tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])
                    filesmeta = scraper._extract_files_meta(tr)
                    for i in filesmeta:
                        rowdict = self.post_search_hook(i, matchdict)
                        rowdict['SatelliteNumber'] = satno
                        rowdict['Level'] = level
                        rowdict['Wavelength'] = wave*u.Angstrom
                        metalist.append(rowdict)

        return QueryResponse(metalist, client=self)

    @classmethod
    def _attrs_module(cls):
        return 'goes', 'sunpy.net.dataretriever.attrs.goes'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        goes_number = [16, 17, 18]
        adict = {attrs.Instrument: [
            ("SUVI", "GOES Solar Ultraviolet Imager.")],
            attrs.goes.SatelliteNumber: [(str(x), f"GOES Satellite Number {x}") for x in goes_number],
            attrs.Source: [('GOES', 'The Geostationary Operational Environmental Satellite Program.')],
            attrs.Physobs: [('flux', 'a measure of the amount of radiation received by an object from a given source.')],
            attrs.Provider: [('NOAA', 'The National Oceanic and Atmospheric Administration.')],
            attrs.Level: [('1b', 'Solar images at six wavelengths with image exposures 10 msec or 1 sec.'),
                          ('2', 'Weighted average of level-1b product files of SUVI.')],
            attrs.Wavelength: [('*')]}
        return adict
