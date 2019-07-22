# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014

from urllib.parse import urljoin

from ..client import GenericClient

__all__ = ['LYRAClient']


class LYRAClient(GenericClient):
    """
    Provides access to the LYRA/Proba2 data `archive <http://proba2.oma.be/lyra/data/bsd/>`__
    hosted by the `PROBA2 Science Center <http://proba2.oma.be>`__.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument('LYRA'))  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the LYRAClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str6     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2016-01-01 00:00:00 2016-01-02 00:00:00 Proba2       lyra        nan
    2016-01-01 00:00:00 2016-01-02 00:00:00 Proba2       lyra        nan
    <BLANKLINE>
    <BLANKLINE>

    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns list of URLS corresponding to value of input timerange.

        Parameters
        ----------
        timerange: sunpy.time.TimeRange
            time range for which data is to be downloaded.

        Returns
        -------
        urls : list
            list of URLs corresponding to the requested time range
        """
        days = timerange.get_dates()
        urls = []
        for day in days:
            urls.append(self._get_url_for_date(day, **kwargs))
        return urls

    def _get_url_for_date(self, date, **kwargs):
        """
        Return URL for corresponding date.

        Parameters
        ----------
        date : `astropy.time.Time`, `~datetime.datetime`, `~datetime.date`

        Returns
        -------
        str
            The URL for the corresponding date.
        """

        filename = "lyra_{}-000000_lev{:d}_std.fits".format(
            date.strftime('%Y%m%d'), kwargs.get('level', 2))
        base_url = "http://proba2.oma.be/lyra/data/bsd/"
        url_path = urljoin(date.strftime('%Y/%m/%d/'), filename)

        return urljoin(base_url, url_path)

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'Proba2'
        self.map_['instrument'] = 'lyra'
        self.map_['physobs'] = 'irradiance'
        self.map_['provider'] = 'esa'

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
        chkattr =  ['Time', 'Instrument', 'Level']
        chklist =  [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'lyra':
                return all(chklist)
        return False
