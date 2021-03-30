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
           Start Time               End Time        Instrument ... Source Provider
    ----------------------- ----------------------- ---------- ... ------ --------
    2016-01-01 00:00:00.000 2016-01-01 23:59:59.999       SOON ...   SWPC     NOAA
    2016-01-02 00:00:00.000 2016-01-02 23:59:59.999       SOON ...   SWPC     NOAA
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
        # Validate time range srs generated daily since 1996
        cur_year = Time.now().datetime.year
        req_start_year = timerange.start.datetime.year
        req_end_year = timerange.end.datetime.year

        # Return early if both start and end are less than or greater than limits
        if req_start_year <= req_end_year < self.MIN_YEAR \
                or req_end_year >= req_start_year > cur_year:
            return result

        # No point in searching below the min or above max years
        start_year = max(req_start_year, self.MIN_YEAR)
        end_year = min(req_end_year, cur_year)

        # Search for tarballs for all years in the query
        tarball_timerange = TimeRange(f'{start_year}-01-01', f'{end_year}-12-31 23:59:59.999')
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
        min_file_date = (datetime(max_tarball_year, 12, 31, 23, 59, 59) if max_tarball_year else
                         datetime(start_year, 1, 1, 0, 0, 0))
        max_file_date = min(timerange.end.datetime, Time.now().datetime)
        if min_file_date < max_file_date:
            file_timerange = TimeRange(f'{min_file_year}-01-01', max_file_date)
            srsfile_scraper = Scraper(self.BASE_URL + '%Y/SRS/%Y%m%dSRS.txt')
            srsfiles = srsfile_scraper.filelist(file_timerange)
            for srs_url in srsfiles:
                date = srsfile_scraper._extractDateURL(srs_url)
                srs_urls[(date.datetime.year, date.datetime.month, date.datetime.day)] = srs_url
                log.debug('SRS file found for year %d', date)

        # Now iterate over all days and if the day is in a year we have a tarball for or a day there
        # is a individual srs file add to the result with corresponding extdict
        for day in timerange.get_dates():
            day_ymd = (int(day.strftime('%Y')), int(day.strftime('%m')), int(day.strftime('%d')))
            extdict = {'year': day_ymd[0], 'month': day_ymd[1], 'day': day_ymd[2]}
            if self.MIN_YEAR <= day_ymd[0] <= cur_year:
                if day_ymd[0] in tarball_urls.keys():
                    result.append((extdict, tarball_urls[day_ymd[0]]))
                elif day_ymd in srs_urls.keys():
                    result.append((extdict, srs_urls[day_ymd]))

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
        urls = [qrblock['url'] for qrblock in qres]
        filenames = []
        local_filenames = []
        for url, qre in zip(urls, qres):
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

        downloader = Downloader(max_conn=2)
        for aurl, fname in zip(urls, paths):
            downloader.enqueue_file(aurl, filename=fname)

        paths = downloader.download()

        outfiles = []
        for fname, srs_filename in zip(local_paths, local_filenames):
            name = fname.name
            past_year = False
            for fname2 in paths:
                fname2 = pathlib.Path(fname2)
                if fname2.name.endswith('.txt'):
                    continue
                year = fname2.name.split('_SRS')[0]
                if year in name:
                    with tarfile.open(fname2) as open_tar:
                        filepath = fname.parent
                        try:
                            member = open_tar.getmember('SRS/' + srs_filename)
                        except KeyError:
                            # Some tars have a {year}_SRS when extracted, 2010 being one example
                            member = open_tar.getmember(f'{year}_SRS/' + srs_filename)
                        member.name = name
                        open_tar.extract(member, path=filepath)
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
