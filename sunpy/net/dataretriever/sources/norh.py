#  Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#  This Module was developed under funding provided by
#  Google Summer of Code 2014

from sunpy.extern.six.moves.urllib.parse import urljoin

from ..client import GenericClient

__all__ = ['NoRHClient']


class NoRHClient(GenericClient):

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns list of URLS corresponding to value of input timerange.

        Parameters
        ----------
        timerange: `sunpy.time.TimeRange`
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
        date : Python datetime object

        Returns
        -------
        string
            The URL for the corresponding date.
        """

        # Hack to get around Python 2.x not backporting PEP 3102.
        wavelength = kwargs.pop('wavelength', None)

        # default urllib password anonymous@ is not accepted by the NoRH FTP
        # server. include an accepted password in base url
        baseurl = 'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/'

        # date is a datetime.date object
        if wavelength == '34':
            final_url = urljoin(baseurl,
                                date.strftime('%Y/%m/' + 'tcz' + '%y%m%d'))
        else:
            final_url = urljoin(baseurl,
                                date.strftime('%Y/%m/' + 'tca' + '%y%m%d'))

        return final_url

    def _makeimap(self):
        """
        Helper Function used to hold information about source.
        """
        self.map_['source'] = 'NAOJ'
        self.map_['provider'] = 'NRO'
        self.map_['instrument'] = 'RadioHelioGraph'
        self.map_['phyobs'] = ''

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
            if x.__class__.__name__ == 'Instrument' and x.value == 'norh':
                return all(chklist)
        return False
