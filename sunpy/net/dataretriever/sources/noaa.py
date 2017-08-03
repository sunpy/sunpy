# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014


from ..client import GenericClient
import datetime
import os
import tarfile
from functools import partial
from collections import OrderedDict

import sunpy
from sunpy.util import replacement_filename, deprecated
from sunpy.net.dataretriever.client import simple_path

from sunpy.net.download import Downloader, Results


__all__ = ['NOAAIndicesClient', 'NOAAPredictClient', 'SRSClient']


class NOAAIndicesClient(GenericClient):

    @staticmethod
    def _get_default_uri():
        """Return the url to download indices"""
        return ["ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt"]

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Helper function:
        """
        return NOAAIndicesClient._get_default_uri()

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'sdic'
        self.map_['instrument'] = 'noaa-indices'
        self.map_['physobs'] = 'sunspot number'
        self.map_['provider'] = 'swpc'

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
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'noaa-indices':
                return all(chklist)
        return False


class NOAAPredictClient(GenericClient):

    @staticmethod
    def _get_default_uri():
        """Return the url to download indices"""
        return ["http://services.swpc.noaa.gov/text/predicted-sunspot-radio-flux.txt"]

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Helper function:
        """
        return NOAAPredictClient._get_default_uri()

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ises'
        self.map_['instrument'] = 'noaa-predict'
        self.map_['physobs'] = 'sunspot number'
        self.map_['provider'] = 'swpc'

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
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'noaa-predict':
                return all(chklist)
        return False


class SRSClient(GenericClient):

    @staticmethod
    def _get_default_uri():
        today = datetime.datetime.utcnow()
        return [('ftp://ftp.swpc.noaa.gov/pub/warehouse/',
                 '{date:%Y}/SRS/{date:%Y%m%d}SRS.txt').format(date=today)]

    def _get_url_for_timerange(self, timerange, **kwargs):

        if not timerange:
            return SRSClient._get_default_uri()
        result = list()
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/warehouse/'
        total_days = (timerange.end - timerange.start).days + 1
        all_dates = timerange.split(total_days)
        today_year = datetime.datetime.utcnow().year
        for day in all_dates:
            if today_year == day.end.year:
                suffix = '{date:%Y}/SRS/{date:%Y%m%d}SRS.txt'
            else:
                suffix = '{date:%Y}/{date:%Y}_SRS.tar.gz'
            url = base_url + suffix.format(date=day.end)
            result.append(url)

        return result

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

        urls = [qrblock.url for qrblock in qres]

        filenames = []
        local_filenames = []

        for i, [url, qre] in enumerate(zip(urls, qres)):
            name = url.split('/')[-1]

            # temporary fix !!! coz All QRBs have same start_time values
            day = qre.time.start.date() + datetime.timedelta(days=i)

            if name not in filenames:
                filenames.append(name)

            if name.endswith('.gz'):
                local_filenames.append('{date:%Y%m%d}SRS.txt'.format(date=day))
            else:
                local_filenames.append(name)

        # Files to be actually downloaded
        paths = self._get_full_filenames(qres, filenames, path)

        # Those files that will be present after get returns
        local_paths = self._get_full_filenames(qres, local_filenames, path)

        res = Results(lambda x: None, 0, lambda map_: self._link(map_))

        # remove duplicate urls. This will make paths and urls to have same number of elements.
        # OrderedDict is required to maintain ordering because it will be zipped with paths later

        urls = list(OrderedDict.fromkeys(urls))

        dobj = Downloader(max_conn=len(urls), max_total=len(urls))

        # We cast to list here in list(zip... to force execution of
        # res.require([x]) at the start of the loop.
        for aurl, ncall, fname in list(zip(urls, map(lambda x: res.require([x]),
                                                     urls), paths)):
            dobj.download(aurl, fname, ncall, error_callback)

        res.wait()

        res2 = Results(lambda x: None, 0)

        for fname, srs_filename in zip(local_paths, local_filenames):

            fname = fname.args[0]
            name = fname.split('/')[-1]

            past_year = False
            for i, fname2 in enumerate(paths):

                fname2 = fname2.args[0]

                if fname2.endswith('.txt'):
                    continue

                year = fname2.split('/')[-1]
                year = year.split('_SRS')[0]

                if year in name:
                    TarFile = tarfile.open(fname2)
                    filepath = fname.rpartition('/')[0]
                    member = TarFile.getmember('SRS/' + srs_filename)
                    member.name = name
                    TarFile.extract(member, path=filepath)
                    TarFile.close()

                    callback = res2.require([fname])
                    callback({'path': fname})

                    past_year = True
                    break

            if past_year is False:
                callback = res2.require([fname])
                callback({'path': fname})

        return res2

    @deprecated('0.8', alternative='NOAAPredictClient.fetch')
    def get(self, qres, path=None, error_callback=None, **kwargs):
        """
        See `~sunpy.net.dataretriever.sources.noaa.NOAAPredictClient.fetch`
        """
        return self.fetch(qres, path=path, error_callback=error_callback, **kwargs)

    def _makeimap(self):
        self.map_['source'] = 'swpc'
        self.map_['instrument'] = 'SOON'
        self.map_['physobs'] = 'SRS'
        self.map_['source'] = 'NOAA/USAF'

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
        chkattr = ["Time", "Instrument"]
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == "Instrument" and\
               str(x.value).lower() in ["soon", "srs_table"]:
                return True
        return False
