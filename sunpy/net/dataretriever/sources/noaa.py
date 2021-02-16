# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014
import pathlib
import tarfile
from datetime import datetime
from collections import OrderedDict

from astropy.time import Time

from sunpy import log
from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.time import TimeRange
from sunpy.util.parfive_helpers import Downloader

__all__ = ['NOAAIndicesClient', 'NOAAPredictClient', 'SRSClient']

from sunpy.util.scraper import Scraper


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

    Uses the `ftp archive <ftp://ftp.swpc.noaa.gov/pub/warehouse/>`__.

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
         Start Time           End Time      Instrument Physobs Source Provider
    ------------------- ------------------- ---------- ------- ------ --------
    2016-01-01 00:00:00 2016-01-01 23:59:59       SOON     SRS   SWPC     NOAA
    2016-01-02 00:00:00 2016-01-02 23:59:59       SOON     SRS   SWPC     NOAA
    <BLANKLINE>
    <BLANKLINE>

    """
    BASE_URL = 'ftp://ftp.swpc.noaa.gov/pub/warehouse/'
    MIN_YEAR = 1996

    def _get_url_for_timerange(self, timerange):
        """
        Returns a list of urls corresponding to a given time-range.
        """
        result = list()

        # Validate time range and return early if out side min/max
        cur_year = Time.now().datetime.year
        req_start_year = timerange.start.datetime.year
        req_end_year = timerange.end.datetime.year
        start_year = req_start_year if self.MIN_YEAR < req_start_year <= cur_year else None
        end_year = req_end_year if self.MIN_YEAR < req_end_year <= cur_year else None
        if start_year is None and end_year is None:
            return result

        # Search for tarballs for all years in the query
        if end_year == start_year and end_year != cur_year:
            end_year += 1
        tarball_timerange = TimeRange(f'{start_year}-01-01', f'{end_year}-01-01')
        tarball_urls = dict()
        tarball_scraper = Scraper(self.BASE_URL + '%Y/%Y_SRS.tar.gz')
        tarballs = tarball_scraper.filelist(tarball_timerange)
        max_tarball_year = None
        for tb_url in tarballs:
            date = tarball_scraper._extractDateURL(tb_url)
            year = date.to_datetime().year
            max_tarball_year = year
            tarball_urls[year] = tb_url
            log.debug('SRS tarball found for year %d', year)

        # Create a new time range for the times not covered by tarballs, have to assume tarballs
        # cover a year, and look for individual srs file for this time range.
        srs_urls = dict()
        min_file_year = max_tarball_year if max_tarball_year else start_year
        min_file_date = datetime(min_file_year, 1, 1, 0, 0, 0)
        if min_file_date <= timerange.end.datetime:
            file_timerange = TimeRange(f'{min_file_year}-01-01', timerange.end)
            srsfile_scraper = Scraper(self.BASE_URL + '%Y/SRS/%Y%m%dSRS.txt')
            srsfiles = srsfile_scraper.filelist(file_timerange)
            for srs_url in srsfiles:
                date = srsfile_scraper._extractDateURL(srs_url)
                srs_urls[(date.datetime.year, date.datetime.month, date.datetime.day)] = srs_url
                log.debug('SRS file found for year %d', date)

        # Now iterate over all days and if the day is in a year we have a tarball for or a day there
        # is a individual srs file add to the result with corresponding extdict
        for day in timerange.get_dates():
            day_ymd = (day.datetime.year, day.datetime.month, day.datetime.day)
            extdict = {name: getattr(day.datetime, name) for name in ['year', 'month', 'day']}
            if self.MIN_YEAR < day_ymd[0] <= cur_year:
                if day_ymd[0] in tarball_urls.keys():
                    result.append((extdict, tarball_urls[day_ymd[0]]))
                elif day_ymd in srs_urls.keys():
                    result.append((extdict, srs_urls[day_ymd]))

        return result

    # def post_search_hook(self, exdict, matchdict):
    #     # update the extracted metadata to include the queried times rather
    #     # than those scraped from the downloaded zip (which includes full year data).
    #     rowdict = super().post_search_hook(exdict, matchdict)
    #     rowdict["Start Time"] = matchdict["Start Time"]
    #     rowdict["End Time"] = matchdict["End Time"]
    #     rowdict["Start Time"].format = 'iso'
    #     rowdict["End Time"].format = 'iso'
    #     return rowdict

    def search(self, *args, **kwargs):
        extractor1 = '{}/warehouse/{:4d}/SRS/{year:4d}{month:2d}{day:2d}SRS.txt'
        extractor2 = '{}/warehouse/{year:4d}/{}'
        matchdict = self._get_match_dict(*args, **kwargs)
        timerange = TimeRange(matchdict['Start Time'], matchdict['End Time'])
        metalist = []
        for extdict, url in self._get_url_for_timerange(timerange):
            extdict['Url'] = url
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

        urls = [qrblock['Url'] for qrblock in qres]

        filenames = []
        local_filenames = []

        for i, [url, qre] in enumerate(zip(urls, qres)):
            name = url.split('/')[-1]

            day = qre['Start Time']

            if name not in filenames:
                filenames.append(name)

            if name.endswith('.gz'):
                local_filenames.append('{}SRS.txt'.format(day.strftime('%Y%m%d')))
            else:
                local_filenames.append(name)

        if path is not None:
            path = pathlib.Path(path)
        # Files to be actually downloaded
        paths = self._get_full_filenames(qres, filenames, path)

        # Those files that will be present after get returns
        local_paths = self._get_full_filenames(qres, local_filenames, path)

        # remove duplicate urls. This will make paths and urls to have same number of elements.
        # OrderedDict is required to maintain ordering because it will be zipped with paths later
        urls = list(OrderedDict.fromkeys(urls))

        dobj = Downloader(max_conn=5)

        for aurl, fname in zip(urls, paths):
            dobj.enqueue_file(aurl, filename=fname)

        paths = dobj.download()

        outfiles = []
        for fname, srs_filename in zip(local_paths, local_filenames):

            name = fname.name

            past_year = False
            for i, fname2 in enumerate(paths):
                fname2 = pathlib.Path(fname2)

                if fname2.name.endswith('.txt'):
                    continue

                year = fname2.name.split('_SRS')[0]

                if year in name:
                    TarFile = tarfile.open(fname2)
                    filepath = fname.parent
                    member = TarFile.getmember('SRS/' + srs_filename)
                    member.name = name
                    TarFile.extract(member, path=filepath)
                    TarFile.close()

                    outfiles.append(fname)

                    past_year = True
                    break

            if past_year is False:
                outfiles.append(fname)

        paths.data = list(map(str, outfiles))
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
