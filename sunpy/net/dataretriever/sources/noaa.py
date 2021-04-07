# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014
import pathlib

from astropy.time import Time

from sunpy import log
from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.time import TimeRange
from sunpy.util.parfive_helpers import Downloader
from sunpy.util.scraper import Scraper

__all__ = ['NOAAIndicesClient', 'NOAAPredictClient', 'SRSClient']


class NOAAIndicesClient(GenericClient):
    """
    Provides access to the NOAA solar cycle indices.

    Uses the `SWPC NOAA archive <https://services.swpc.noaa.gov/json/solar-cycle/>`__.
    This is a fixed dataset so the result is independent of the time range.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.noaa_indices)  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the NOAAIndicesClient:
     Instrument     Physobs     Source Provider
    ------------ -------------- ------ --------
    NOAA-Indices sunspot number   SIDC     SWPC
    <BLANKLINE>
    <BLANKLINE>

    """
    required = {a.Instrument}

    def search(self, *args, **kwargs):
        rowdict = self._get_match_dict(*args, **kwargs)
        for key in rowdict:
            if isinstance(rowdict[key], list):
                # uses first value among the list of possible values corresponding to an Attr
                # returned by `get_match_dict()` to be shown in query response table.
                rowdict[key] = rowdict[key][0]
        rowdict['url'] = 'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json'
        rowdict['Instrument'] = 'NOAA-Indices'
        # These results are not dependant on time, but we allow time as a
        # parameter for easy searching, so remove time from the results table
        # injected by GenericClient.
        rowdict.pop('Start Time', None)
        rowdict.pop('End Time', None)
        return QueryResponse([rowdict], client=self)

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [
            ('NOAA-Indices', 'Recent Solar Indices of Observed Monthly Mean Values')],
            attrs.Physobs: [('sunspot number', 'Sunspot Number.')],
            attrs.Source: [('SIDC', 'The Solar Influence Data Analysis Center')],
            attrs.Provider: [('SWPC', 'The Space Weather Prediction Center.')],
            attrs.Time: [('*')]}
        return adict


class NOAAPredictClient(GenericClient):
    """
    Provides access to the NOAA SWPC predicted sunspot Number and 10.7 cm radio flux values.

    Uses the `SWPC NOAA archive <https://services.swpc.noaa.gov/json/solar-cycle/>`__.
    This is a fixed prediction so the result is independent of the time range.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.noaa_predict)  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the NOAAPredictClient:
     Instrument     Physobs     Source Provider
    ------------ -------------- ------ --------
    NOAA-Predict sunspot number   ISES     SWPC
    <BLANKLINE>
    <BLANKLINE>

    """
    required = {a.Instrument}

    def search(self, *args, **kwargs):
        rowdict = self._get_match_dict(*args, **kwargs)
        for key in rowdict:
            if isinstance(rowdict[key], list):
                # uses first value among the list of possible values corresponding to an Attr
                # returned by `get_match_dict()` to be shown in query response table.
                rowdict[key] = rowdict[key][0]
        rowdict['url'] = 'https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json'
        rowdict['Instrument'] = 'NOAA-Predict'
        # These results are not dependant on time, but we allow time as a
        # parameter for easy searching, so remove time from the results table
        # injected by GenericClient.
        rowdict.pop('Start Time', None)
        rowdict.pop('End Time', None)
        return QueryResponse([rowdict], client=self)

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [
            ('NOAA-Predict', 'Predicted Sunspot Number And Radio Flux Values With Expected Ranges.')],
            attrs.Physobs: [('sunspot number', 'Sunspot Number.')],
            attrs.Source: [('ISES', 'The International Space Environmental Services.')],
            attrs.Provider: [('SWPC', 'The Space Weather Prediction Center.')],
            attrs.Time: [('*')]}
        return adict


class SRSClient(GenericClient):
    """
    Provides access to the NOAA SWPC solar region summary data.

    Uses the `ftp archive <https://www.ngdc.noaa.gov/stp/spaceweather.html>`__.

    Notes
    -----
    Data pre-1996 is in free-form text, which cannot be parsed by sunpy, and
    therefore only results from 1996 onwards are returned by this client.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.soon)  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the SRSClient:
           Start Time               End Time        Instrument ... Source Provider
    ----------------------- ----------------------- ---------- ... ------ --------
    2016-01-01 00:00:00.000 2016-01-01 23:59:59.999       SOON ...   SWPC     NOAA
    2016-01-02 00:00:00.000 2016-01-02 23:59:59.999       SOON ...   SWPC     NOAA
    <BLANKLINE>
    <BLANKLINE>

    """
    BASE_URL = 'ftp://ftp.ngdc.noaa.gov/STP/swpc_products/daily_reports/solar_region_summaries/'
    MIN_YEAR = 1996

    def _get_url_for_timerange(self, timerange):
        """
        Returns a list of urls corresponding to a given time-range.
        """
        result = list()
        # Validate time range srs generated daily since 1996
        cur_year = Time.now().datetime.year
        req_start_year = timerange.start.datetime.year
        req_end_year = timerange.end.datetime.year

        # Return early if both start and end are less than or greater than limits
        if req_start_year <= req_end_year < self.MIN_YEAR or req_end_year >= req_start_year > cur_year:
            return result

        # Correct given timerange
        min_file_year = self.MIN_YEAR if self.MIN_YEAR > req_start_year else req_start_year
        min_file_date = Time(f"{min_file_year}-01-01")
        new_start = max(timerange.start.datetime, min_file_date)
        new_end = min(timerange.end.datetime, Time.now().datetime)
        timerange = TimeRange(new_start, new_end)

        srsfile_scraper = Scraper(self.BASE_URL + '%Y/%m/%Y%m%dSRS.txt')
        srsfiles = srsfile_scraper.filelist(timerange)
        srs_urls = dict()
        for srs_url in srsfiles:
            date = srsfile_scraper._extractDateURL(srs_url)
            srs_urls[(date.datetime.year, date.datetime.month, date.datetime.day)] = srs_url
            day_ymd = (int(date.strftime('%Y')), int(date.strftime('%m')), int(date.strftime('%d')))
            extdict = {'year': day_ymd[0], 'month': day_ymd[1], 'day': day_ymd[2]}
            result.append((extdict, srs_urls[day_ymd]))
            log.debug('SRS file found for year %d', date)
        return result

    def search(self, *args, **kwargs):
        matchdict = self._get_match_dict(*args, **kwargs)
        timerange = TimeRange(matchdict['Start Time'], matchdict['End Time'])
        metalist = []
        for extdict, url in self._get_url_for_timerange(timerange):
            extdict['url'] = url
            rowdict = self.post_search_hook(extdict, matchdict)
            metalist.append(rowdict)
        return QueryResponse(metalist, client=self)

    def fetch(self, qres, path=None, error_callback=None, **kwargs):
        """
        Download a set of results.

        Parameters
        ----------
        qres : `~sunpy.net.dataretriever.QueryResponse`
            Results to download.

        Returns
        -------
        Results Object
        """
        if path is not None:
            path = pathlib.Path(path)
        urls = sorted([qrblock['url'] for qrblock in qres])
        filenames = [pathlib.Path(url).name for url in urls]
        # Files to be actually downloaded
        paths = self._get_full_filenames(qres, filenames, path)
        downloader = Downloader(max_conn=2)
        for aurl, fname in zip(urls, paths):
            # Need to change the passive command as the server does not support the aioftp default
            downloader.enqueue_file(aurl, filename=fname, passive_commands=["pasv"])
        paths = downloader.download()
        return paths

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [("SOON", "Solar Region Summary."),
                                    ("SRS-Table", "Solar Region Summary.")],
                 attrs.Physobs: [('SRS', 'Solar Region Summary.')],
                 attrs.Source: [('SWPC', 'The Space Weather Prediction Center.')],
                 attrs.Provider: [('NOAA', 'The National Oceanic and Atmospheric Administration.')]}
        return adict
