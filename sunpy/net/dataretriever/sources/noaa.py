# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014
import pathlib
import tarfile
from collections import OrderedDict

from astropy import units as u
from astropy.time import Time, TimeDelta

from sunpy.extern.parse import parse
from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.util.parfive_helpers import Downloader

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
    2016-01-01 00:00:00 2016-12-31 23:59:59       SOON     SRS   SWPC     NOAA
    2016-01-01 00:00:00 2016-12-31 23:59:59       SOON     SRS   SWPC     NOAA
    <BLANKLINE>
    <BLANKLINE>

    """

    def _get_url_for_timerange(self, timerange):
        """
        Returns a list of urls corresponding to a
        given time-range.
        """
        result = list()
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/warehouse/'
        total_days = int(timerange.days.value) + 1
        all_dates = timerange.split(total_days)
        today_year = int(Time.now().strftime('%Y'))
        for day in all_dates:
            end_year = int(day.end.strftime('%Y'))
            if end_year > today_year or end_year < 1996:
                continue
            elif end_year == today_year:
                suffix = '{}/SRS/{}SRS.txt'.format(
                    end_year, day.end.strftime('%Y%m%d'))
            else:
                suffix = '{}/{}_SRS.tar.gz'.format(
                    end_year, day.end.strftime('%Y'))
            url = base_url + suffix
            result.append(url)
        return result

    def post_search_hook(self, exdict, matchdict):
        # update the extracted metadata to include the queried times rather
        # than those scraped from the downloaded zip (which includes full year data).
        rowdict = super().post_search_hook(exdict, matchdict)
        rowdict["Time"] = matchdict["Time"]
        return rowdict

    def search(self, *args, **kwargs):
        extractor1 = '{}/warehouse/{:4d}/SRS/{year:4d}{month:2d}{day:2d}SRS.txt'
        extractor2 = '{}/warehouse/{year:4d}/{}'
        matchdict = self._get_match_dict(*args, **kwargs)
        timerange = matchdict['Time']
        metalist = []
        for url in self._get_url_for_timerange(timerange):
            exdict1 = parse(extractor1, url)
            exdict2 = parse(extractor2, url)
            exdict = (exdict2 if exdict1 is None else exdict1).named
            exdict['url'] = url
            rowdict = self.post_search_hook(exdict, matchdict)
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

        for i, [url, qre] in enumerate(zip(urls, qres)):
            name = url.split('/')[-1]

            day = Time(qre['Time'].start.strftime('%Y-%m-%d')) + TimeDelta(i*u.day)

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
