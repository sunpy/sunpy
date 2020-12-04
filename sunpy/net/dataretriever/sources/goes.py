# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014

from datetime import datetime
from collections import OrderedDict

import astropy.units as u

from sunpy import config
from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.time import TimeRange
from sunpy.util.scraper import Scraper, get_timerange_from_exdict

TIME_FORMAT = config.get("general", "time_format")

__all__ = ["XRSClient", "SUVIClient"]


class XRSClient(GenericClient):
    """
    Provides access to the GOES XRS fits files archive.

    Searches for GOES XRS data both on NASA servers prior to re-processed
    GOES 13, 14 and 15 and on the NOAA archive for > GOES 13.
    For satellite numbers > 13 the XRSClient searches the NOAA archive, and
    returns the re-processed science-quality data for GOES 13, 14 and 15, and
    also the new GOES-R series 16 and 17.

    Note - the new science quality data have scaling factors removed for 13, 14 and 15
    and they are not added to GOES 16 AND 17. This means the peak flux will be different to
    the older version of the data, such as those collected from the NASA servers.

    See the following readmes about the data

    * Reprocessed 13, 14, 15 :
        https://satdat.ngdc.noaa.gov/sem/goes/data/science/xrs/GOES_13-15_XRS_Science-Quality_Data_Readme.pdf

    * GOES-R 16, 17 :
         https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l1b/docs/GOES-R_XRS_L1b_Science-Quality_Data_Readme.pdf

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
         Start Time           End Time      Instrument ... Source Provider
    ------------------- ------------------- ---------- ... ------ --------
    2016-01-01 00:00:00 2016-01-01 23:59:59        XRS ...   GOES     NOAA
    2016-01-02 00:00:00 2016-01-02 23:59:59        XRS ...   GOES     NOAA
    2016-01-01 00:00:00 2016-01-01 23:59:59        XRS ...   GOES     NOAA
    2016-01-02 00:00:00 2016-01-02 23:59:59        XRS ...   GOES     NOAA
    <BLANKLINE>
    <BLANKLINE>
    """
    # GOES XRS data from NASA servers upto GOES 15. The reprocessed 13, 14, 15 data should be taken from NOAA server.
    baseurl_old = r'https://umbra.nascom.nasa.gov/goes/fits/%Y/go(\d){2}(\d){6,8}\.fits'
    pattern_old = '{}/fits/{year:4d}/go{SatelliteNumber:02d}{}{month:2d}{day:2d}.fits'

    # GOES XRS 13, 14, 15 from NOAA (re-processed data)
    baseurl_new = (r"https://satdat.ngdc.noaa.gov/sem/goes/data/science/xrs/"
                   r"goes{SatelliteNumber}/gxrs-l2-irrad_science/%Y/%m/sci_gxrs-l2-irrad_g{SatelliteNumber}_d%Y%m%d_.*\.nc")
    pattern_new = ("{}/goes{SatelliteNumber:02d}/gxrs-l2-irrad_science/{year:4d}/"
                   "{month:2d}/sci_gxrs-l2-irrad_g{SatelliteNumber:02d}_d{year:4d}{month:2d}{day:2d}_{}.nc")

    # GOES XRS data for GOES-R Series - 16, 17
    baseurl_r = (r"https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes{SatelliteNumber}"
                 r"/l2/data/xrsf-l2-flx1s_science/%Y/%m/sci_xrsf-l2-flx1s_g{SatelliteNumber}_d%Y%m%d_.*\.nc")
    pattern_r = ("{}/goes/goes{SatelliteNumber:02d}/l2/data/xrsf-l2-flx1s_science/{year:4d}/"
                 "{month:2d}/sci_xrsf-l2-flx1s_g{SatelliteNumber:02d}_d{year:4d}{month:2d}{day:2d}_{}.nc")

    def post_search_hook(self, i, matchdict):
        tr = get_timerange_from_exdict(i)
        rowdict = OrderedDict()
        rowdict["Time"] = TimeRange(tr.start, tr.end)
        rowdict["Start Time"] = tr.start.strftime(TIME_FORMAT)
        rowdict["End Time"] = tr.end.strftime(TIME_FORMAT)
        rowdict["Instrument"] = matchdict["Instrument"][0].upper()
        rowdict["SatelliteNumber"] = i["SatelliteNumber"]
        rowdict["Physobs"] = matchdict["Physobs"][0]
        rowdict["url"] = i["url"]
        rowdict["Source"] = matchdict["Source"][0]
        if i["url"].endswith(".fits"):
            rowdict["Provider"] = matchdict["Provider"][0]
        else:
            rowdict["Provider"] = matchdict["Provider"][1]

        return rowdict

    def search(self, *args, **kwargs):

        matchdict = self._get_match_dict(*args, **kwargs)

        # this is for the case when the timerange overlaps with the provider change.
        if matchdict["Time"].start < "2009-09-01" and matchdict["Time"].end >= "2009-09-01":
            matchdict_before, matchdict_after = matchdict.copy(), matchdict.copy()
            matchdict_before["Time"] = TimeRange(matchdict_before["Time"].start, '2009-08-31')
            matchdict_after["Time"] = TimeRange('2009-09-01', matchdict_after["Time"].end)
            metalist_before = self._get_metalist(matchdict_before)
            metalist_after = self._get_metalist(matchdict_after)
            metalist = metalist_before + metalist_after
        else:
            metalist = self._get_metalist(matchdict)

        return QueryResponse(metalist, client=self)

    def _get_metalist_fn(self, matchdict, baseurl, pattern):
        """
        Function to help get list of OrderedDicts.
        """
        metalist = []
        scraper = Scraper(baseurl, regex=True)
        filemeta = scraper._extract_files_meta(
            matchdict["Time"], extractor=pattern, matcher=matchdict)
        for i in filemeta:
            rowdict = self.post_search_hook(i, matchdict)
            metalist.append(rowdict)
        return metalist

    def _get_metalist(self, matchdict):
        """
        Function to get the list of OrderDicts.
        This makes it easier for when searching for overlapping providers.
        """
        metalist = []
        # the data before the re-processed GOES 13, 14, 15 data.
        if (matchdict["Time"].end < "2009-09-01") or (matchdict["Time"].end >= "2009-09-01" and matchdict["Provider"] == ["sdac"]):
            metalist += self._get_metalist_fn(matchdict, self.baseurl_old, self.pattern_old)

        # new data from NOAA.
        else:
            if matchdict["Time"].end >= "2017-02-07":
                for sat in [16, 17]:
                    metalist += self._get_metalist_fn(matchdict,
                                                      self.baseurl_r.format(SatelliteNumber=sat), self.pattern_r)

            if matchdict["Time"].end <= "2020-03-04":
                for sat in [13, 14, 15]:
                    metalist += self._get_metalist_fn(matchdict,
                                                      self.baseurl_new.format(SatelliteNumber=sat), self.pattern_new)

        return metalist

    @classmethod
    def _attrs_module(cls):
        return 'goes', 'sunpy.net.dataretriever.attrs.goes'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        goes_number = [2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
        adict = {attrs.Instrument: [
            ("GOES", "The Geostationary Operational Environmental Satellite Program."),
            ("XRS", "GOES X-ray Sensor")],
            attrs.Physobs: [('irradiance', 'the flux of radiant energy per unit area.')],
            attrs.Source: [("GOES", "The Geostationary Operational Environmental Satellite Program.")],
            attrs.Provider: [('SDAC', 'The Solar Data Analysis Center.'),
                             ('NOAA', 'The National Oceanic and Atmospheric Administration.')],
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

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> import astropy.units as u
    >>> results = Fido.search(a.Time("2020/7/10", "2020/7/10 00:10"), a.Instrument('suvi'),a.Level.two,
    ...                       a.goes.SatelliteNumber(16), a.Wavelength(304*u.Angstrom))  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the SUVIClient:
         Start Time           End Time      Instrument ... Level   Wavelength
    ------------------- ------------------- ---------- ... ----- --------------
    2020-07-10 00:00:00 2020-07-10 00:04:00       SUVI ...     2 304.0 Angstrom
    2020-07-10 00:04:00 2020-07-10 00:08:00       SUVI ...     2 304.0 Angstrom
    2020-07-10 00:08:00 2020-07-10 00:12:00       SUVI ...     2 304.0 Angstrom
    <BLANKLINE>
    <BLANKLINE>

    """
    baseurl1b = (r'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes'
                 r'{SatelliteNumber}/l1b/suvi-l1b-{elem:2}{wave:03}/%Y/%m/%d/OR_SUVI-L1b.*\.fits.gz')
    pattern1b = ('{}/goes/goes{SatelliteNumber:2d}/l{Level:2w}/suvi-l1b-{}{Wavelength:03d}/'
                 '{year:4d}/{month:2d}/{day:2d}/{}_s{:7d}{hour:2d}{minute:2d}{second:2d}'
                 '{:1d}_e{:7d}{ehour:2d}{eminute:2d}{esecond:2d}{:1d}_{}')
    baseurl2 = (r'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes{SatelliteNumber}/'
                r'l2/data/suvi-l2-ci{wave:03}/%Y/%m/%d/dr_suvi-l2-ci{wave:03}_g{SatelliteNumber}_s%Y%m%dT%H%M%SZ_.*\.fits')
    pattern2 = ('{}/goes/goes{SatelliteNumber:2d}/{}/dr_suvi-l{Level}-ci{Wavelength:03d}_g{SatelliteNumber:2d}_s'
                '{year:4d}{month:2d}{day:2d}T{hour:2d}{minute:2d}{second:2d}Z_e'
                '{eyear:4d}{emonth:2d}{eday:2d}T{ehour:2d}{eminute:2d}{esecond:2d}Z_{}')

    def post_search_hook(self, i, matchdict):

        # extracting start times and end times
        start = datetime(i['year'], i['month'], i['day'], i['hour'], i['minute'], i['second'])
        end = datetime(i['year'], i['month'], i['day'], i['ehour'], i['eminute'], i['esecond'])
        timerange = TimeRange(start, end)

        rowdict = OrderedDict()
        rowdict['Time'] = timerange
        rowdict['Start Time'] = start.strftime(TIME_FORMAT)
        rowdict['End Time'] = end.strftime(TIME_FORMAT)
        rowdict['Instrument'] = matchdict['Instrument'][0].upper()
        rowdict['Physobs'] = matchdict['Physobs'][0]
        rowdict['Source'] = matchdict['Source'][0]
        rowdict['Provider'] = matchdict['Provider'][0]
        rowdict['SatelliteNumber'] = i['SatelliteNumber']
        rowdict['Level'] = i['Level']
        rowdict['Wavelength'] = i['Wavelength']*u.Angstrom
        rowdict['url'] = i['url']
        return rowdict

    def search(self, *args, **kwargs):
        supported_waves = [94, 131, 171, 195, 284, 304]*u.Angstrom
        all_waves = []
        matchdict = self._get_match_dict(*args, **kwargs)
        req_wave = matchdict.get('Wavelength', None)
        if req_wave is not None:
            wmin = req_wave.min.to(u.Angstrom, equivalencies=u.spectral())
            wmax = req_wave.max.to(u.Angstrom, equivalencies=u.spectral())
            req_wave = a.Wavelength(wmin, wmax)
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
                        baseurl = self.baseurl1b
                        pattern = self.pattern1b
                    elif str(level) == '2':
                        baseurl = self.baseurl2
                        pattern = self.pattern2
                    else:
                        raise ValueError(f"Level {level} is not supported.")
                    # formatting baseurl using Level, SatelliteNumber and Wavelength
                    urlpattern = baseurl.format(**formdict)

                    scraper = Scraper(urlpattern)
                    filesmeta = scraper._extract_files_meta(matchdict['Time'], extractor=pattern)
                    for i in filesmeta:
                        rowdict = self.post_search_hook(i, matchdict)
                        metalist.append(rowdict)

        return QueryResponse(metalist, client=self)

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
            attrs.Provider: [('NOAA', 'The National Oceanic and Atmospheric Administration.')],
            attrs.Level: [('1b', 'Solar images at six wavelengths with image exposures 10 msec or 1 sec.'),
                          ('2', 'Weighted average of level-1b product files of SUVI.')],
            attrs.Wavelength: [('*')]}
        return adict
