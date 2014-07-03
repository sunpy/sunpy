from sunpy.net.vso.attrs import Time,Instrument,Level
from sunpy.net.unifieddownloader.client import GenericClient
import datetime,urlparse

__all__ = ['Time','Instrument','Level']

class LYRAClient(GenericClient):

        def _get_url_for_timerange(cls, timerange, **kwargs):
            """
            Helper function:
            Input:
            timerange: Time-range over which data is to be downloaded
	    Output: List of urls
            """               
	    days = timerange.get_days()
            urls = []
	    for day in days:
                urls.append(cls._get_url_for_date(day, **kwargs))
            return urls

        def _get_url_for_date(cls, date, **kwargs):
            """Returns a URL to the LYRA data for the specified date"""
            if not isinstance(date, datetime.date):
                raise ValueError("This method requires a date")
	    filename = "lyra_%s000000_lev%d_%s.fits" % (date.strftime('%Y%m%d-'),2, 'std')
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
        def _can_handle_query(cls,*query):
            """Boolean Function:Answers whether client can service the query.
            """
	    chkattr =  ['Time','Instrument','Level']
            chklist =  [x.__class__.__name__ in chkattr for x in query]
            for x in query:
	        if x.__class__.__name__ == 'Instrument' and x.value == 'lyra':
                    return all(chklist)
	    return False
 

