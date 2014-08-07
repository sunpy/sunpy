import datetime
import urlparse

from sunpy.net.unifieddownloader.client import GenericClient

__all__ = ['NoRHClient']
class NoRHClient(GenericClient):

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns list of URLS corresponding to TimeRange.
        Parameters
	----------
        timerange: TimeRange for which data is to be downloaded.
        Returns
	-------
	List of urls
        """
        if not timerange:
            return []

        days = timerange.get_dates()
        urls = []
        for day in days:
            urls.append(self._get_url_for_date(day, **kwargs))
        return urls

    def _get_url_for_date(self, date, **kwargs):
        """Return URL for corresponding date.
	Parameters
	----------
	date : datetime 

        Returns
	-------
	string representing URL
        """
        # Hack to get around Python 2.x not backporting PEP 3102.
        wavelength = kwargs.pop('wavelength', None)

        #default urllib password anonymous@ is not accepted by the NoRH FTP server.
        #include an accepted password in base url
        baseurl = 'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/'

        #date is a datetime.date object
        if wavelength == '34':
            final_url = urlparse.urljoin(baseurl, date.strftime('%Y/%m/' + 'tcz' + '%y%m%d'))
        else:
            final_url = urlparse.urljoin(baseurl, date.strftime('%Y/%m/' + 'tca' + '%y%m%d'))

        return final_url

    def _makeimap(self):
        '''Helper Function:used to hold information about source. '''
        self.map_['source'] = 'NAOJ'
        self.map_['provider'] ='NRO'
        self.map_['instrument'] = 'RadioHelioGraph'
        self.map_['phyobs'] = ''

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Boolean Function:Answers whether client can service the query.
        """
        chkattr =  ['Time', 'Instrument']
        chklist =  [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'norh':
                return all(chklist)
        return False
