# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014
import socket
import posixpath
from datetime import datetime
from urllib.error import URLError
from urllib.request import urlopen, urlretrieve

from dateutil.rrule import rrule, MONTHLY

import astropy.units as u

from sunpy.time import TimeRange, parse_time
from sunpy.instr import rhessi

from ..client import GenericClient

__all__ = ['RHESSIClient']

data_servers = ('https://hesperia.gsfc.nasa.gov/hessidata/',
                'http://hessi.ssl.berkeley.edu/hessidata/',
                'http://soleil.i4ds.ch/hessidata/')

lc_linecolors = ('black', 'pink', 'green', 'blue', 'brown', 'red',
                 'navy', 'orange', 'green')


def get_base_url():
    """
    Find the first mirror which is online
    """
    for server in data_servers:
        try:
            urlopen(server, timeout=1)
            return server
        except (URLError, socket.timeout):
            pass

    raise IOError('Unable to find an online HESSI server from {0}'.format(data_servers))


class RHESSIClient(GenericClient):

    def get_observing_summary_filename(self, time_range):
        """
        Download the RHESSI observing summary data from one of the RHESSI
        servers, parses it, and returns the name of the observing summary files
        relevant for the time range.

        Parameters
        ----------
        time_range : str, TimeRange
            A TimeRange or time range compatible string

        Returns
        -------
        out : list
            Returns the urls of the observation summary file

        Examples
        --------
        >>> from sunpy.net.dataretriever.sources.rhessi import RHESSIClient
        >>> RHESSIClient().get_observing_summary_filename(('2011/04/04', '2011/04/04'))   # doctest: +REMOTE_DATA
        ['https://.../hessidata/metadata/catalog/hsi_obssumm_20110404_042.fits']
        """
        dt = TimeRange(time_range)
        # remove time from dates
        dt = TimeRange(dt.start.strftime('%Y-%m-%d'), dt.end.strftime('%Y-%m-%d'))

        filenames = []

        diff_months = (dt.end.datetime.year - dt.start.datetime.year) * 12 + dt.end.datetime.month - dt.start.datetime.month
        first_month = datetime(dt.start.datetime.year, dt.start.datetime.month, 1)
        month_list = rrule(MONTHLY, dtstart=first_month, count=diff_months+1)

        # need to download and inspect the dbase file to determine the filename
        # for the observing summary data
        # the dbase files are monthly but contain the daily filenames
        for this_month in month_list:
            dbase_file_name, hdrs = self.get_observing_summary_dbase_file(this_month)
            dbase_dat = rhessi.parse_observing_summary_dbase_file(dbase_file_name)
            this_month_obssumm_filenames = dbase_dat.get('filename')
            daily_filenames_dates = [datetime.strptime(d[0:20], 'hsi_obssumm_%Y%m%d') for d in this_month_obssumm_filenames]
            for i, this_date in enumerate(daily_filenames_dates):
                if dt.start <= this_date <= dt.end:
                    filenames.append(posixpath.join(get_base_url(), 'metadata', 'catalog', this_month_obssumm_filenames[i] + 's'))

        return filenames

    @staticmethod
    def get_observing_summary_dbase_file(time):
        """
        Download the RHESSI observing summary database file for the time given.
        One file covers an entire month.  This file lists the name of observing
        summary files for specific times.

        Parameters
        ----------
        time : `str`, datetime

        Returns
        -------
        value : `tuple`
            Return a `tuple` (filename, headers) where filename is the local file
            name under which the object can be found, and headers is
            whatever the info() method of the object returned by urlopen.

        Examples
        --------
        >>> from sunpy.net.dataretriever.sources.rhessi import RHESSIClient
        >>> fname, headers = RHESSIClient.get_observing_summary_dbase_file('2011/04/04')   # doctest: +REMOTE_DATA

        References
        ----------
        | http://hesperia.gsfc.nasa.gov/ssw/hessi/doc/guides/hessi_data_access.htm#Observing Summary Data

        .. note::
            This API is currently limited to providing data from whole days only.

        """
        _time = parse_time(time)

        if _time < parse_time("2002/02/01"):
            raise ValueError("RHESSI summary files are not available before 2002-02-01")

        url = posixpath.join(get_base_url(), 'dbase', _time.strftime("hsi_obssumm_filedb_%Y%m.txt"))
        return urlretrieve(url)

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns a URL to the RHESSI data for the specified date range.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            Date range should be specified using a TimeRange.
        """
        return self.get_observing_summary_filename(timerange)

    def _get_time_for_url(self, urls):
        ts = [datetime.strptime(url.split("hsi_obssumm_")[1].split("_")[0],
                                "%Y%m%d") for url in urls]
        return [TimeRange(t, (1*u.day-1*u.ms)) for t in ts]

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'rhessi'
        self.map_['instrument'] = 'rhessi'
        self.map_['physobs'] = 'irradiance'
        self.map_['provider'] = 'nasa'

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
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'rhessi':
                return all(chklist)
        return False
