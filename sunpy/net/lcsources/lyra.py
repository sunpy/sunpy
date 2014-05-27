from sunpy.net.lcsources.client import GenericClient
import datetime,urlparse
class LYRAClient(GenericClient):

        def _get_url_for_timerange(cls, timerange, **kwargs):
            
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

        def _makeimap(self,*args,**kwargs):
	    '''map_:Dict'''
            GenericClient.makeargs(self,args,kwargs)
	    self.map_['source'] = 'Proba2'
	    self.map_['instrument'] = 'lyra'
	    self.map_['phyobs'] = 'irradiance'
	    self.map_['provider'] = 'esa'


