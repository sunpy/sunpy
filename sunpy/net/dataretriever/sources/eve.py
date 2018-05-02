# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding by
# Google Summer of Code 2014

from __future__ import absolute_import, division, print_function
import datetime

from sunpy.time import TimeRange
from sunpy.util.scraper import Scraper
from sunpy.extern.six.moves.urllib.parse import urljoin

from ..client import GenericClient


__all__ = ['EVEClient']

BASEURL = ('http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/'
           'L0CS/SpWx/%Y/%Y%m%d_EVE_L0CS_DIODES_1m.txt')


class EVEClient(GenericClient):
    """
    This EVEClient is for the Level 0C data
    from http://lasp.colorado.edu/home/eve/data/data-access/.

    To use this client you must request Level 0 data.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument('EVE'), a.Level(0))  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the EVEClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str3     str3       str3
    ------------------- ------------------- ------ ---------- ----------
    2016-01-01 00:00:00 2016-01-02 00:00:00    SDO        eve        nan
    2016-01-02 00:00:00 2016-01-03 00:00:00    SDO        eve        nan
    <BLANKLINE>
    <BLANKLINE>

    """

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Return list of URLS corresponding to value of input timerange.

        Parameters
        ----------
        timerange: `sunpy.time.TimeRange`
            time range for which data is to be downloaded.

        Returns
        -------
        urls : list
            list of URLs corresponding to the requested time range

        """
        # If start of time range is before 00:00, converted to such, so
        # files of the requested time ranger are included.
        # This is done because the archive contains daily files.
        if timerange.start.time() != datetime.time(0, 0):
            timerange = TimeRange('{:%Y-%m-%d}'.format(timerange.start), timerange.end)
        eve = Scraper(BASEURL)
        return eve.filelist(timerange)

    def _get_time_for_url(self, urls):
        eve = Scraper(BASEURL)
        times = list()
        for url in urls:
            t0 = eve._extractDateURL(url)
            # hard coded full day as that's the normal.
            times.append(TimeRange(t0, t0 + datetime.timedelta(days=1)))
        return times

    def _makeimap(self):
        """
        Helper Function: used to hold information about source.
        """
        self.map_['source'] = 'SDO'
        self.map_['provider'] = 'LASP'
        self.map_['instrument'] = 'eve'
        self.map_['physobs'] = 'irradiance'

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
        chk_var = 0
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'eve':
                chk_var += 1

            elif x.__class__.__name__ == 'Level' and x.value == 0:
                chk_var += 1

        if chk_var == 2:
            return True
        return False
