# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014
import tarfile
import pathlib
from collections import OrderedDict

from parfive import Downloader

import astropy.units as u
from astropy.time import Time, TimeDelta

from ..client import GenericClient

__all__ = ['NOAAIndicesClient', 'NOAAPredictClient', 'SRSClient']


class NOAAIndicesClient(GenericClient):
    @staticmethod
    def _get_url_for_timerange(timerange, **kwargs):
        """
        Helper function:
        """
        return ["ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt"]

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
    def _get_url_for_timerange(self, timerange, **kwargs):
        result = list()
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/warehouse/'
        total_days = int(timerange.days.value) + 1
        all_dates = timerange.split(total_days)
        today_year = Time.now().strftime('%Y')
        for day in all_dates:
            if today_year == day.end.strftime('%Y'):
                suffix = '{0}/SRS/{1}SRS.txt'.format(
                    day.end.strftime('%Y'), day.end.strftime('%Y%m%d'))
            else:
                suffix = '{0}/{1}_SRS.tar.gz'.format(
                    day.end.strftime('%Y'), day.end.strftime('%Y'))
            url = base_url + suffix
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
            day = Time(qre.time.start.strftime('%Y-%m-%d')) + TimeDelta(i*u.day)

            if name not in filenames:
                filenames.append(name)

            if name.endswith('.gz'):
                local_filenames.append('{}SRS.txt'.format(day.strftime('%Y%m%d')))
            else:
                local_filenames.append(name)

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
        for x in query:
            if (x.__class__.__name__ == "Instrument" and
                str(x.value).lower() in ["soon", "srs_table"]):
                return True
        return False
