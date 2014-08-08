#Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#This module was developed under funding provided by 
#Google Summer of Code 2014

import datetime
import urlparse

from sunpy.net.unifieddownloader.client import GenericClient

__all__ = ['LYRAClient']

class LYRAClient(GenericClient):
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
        """
        Return URL for corresponding date.
	Parameters
	----------
	date : datetime 

        Returns
	-------
	string representing URL
	"""
        if not isinstance(date, datetime.date):
            raise ValueError("This method requires a date")
        filename = "lyra_%s000000_lev%d_%s.fits" % (date.strftime('%Y%m%d-'), 2, 'std')
        base_url = "http://proba2.oma.be/lyra/data/bsd/"
        url_path = urlparse.urljoin(date.strftime('%Y/%m/%d/'), filename)
        return urlparse.urljoin(base_url, url_path)

    def _makeimap(self):
        '''Helper Function:used to hold information about source. '''
        self.map_['source'] = 'Proba2'
        self.map_['instrument'] = 'lyra'
        self.map_['phyobs'] = 'irradiance'
        self.map_['provider'] = 'esa'

    @classmethod
    def _can_handle_query(cls, *query):
        """Boolean Function:Answers whether client can service the query.
        """
        chkattr =  ['Time', 'Instrument', 'Level']
        chklist =  [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'lyra':
                return all(chklist)
        return False
