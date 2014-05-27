import urlparse
from sunpy.net.vso.attrs import Time
from sunpy.net.lcsources.client import GenericClient

class EVEClient(GenericClient):

    def _get_url_for_timerange(cls,timerange,**kwargs):
         if not timerange:
              return []
         days = timerange.get_days()
	 urls = []
	 for day in days:
	     urls.append(cls._get_url_for_date(day,**kwargs))
	 return urls

    def _get_url_for_date(cls,date,**kwargs):
        
    #	if date < datetime.date(2010,1,1):
    #        raise ERROR (Define error class showing data not available for this date
        base_url = 'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/'
        return urlparse.urljoin(base_url, date.strftime('%Y/%Y%m%d') + '_EVE_L0CS_DIODES_1m.txt')
    
    def makeimap(self,*args,**kwargs):
       '''map_:Dict'''
       GenericClient.makeargs(self,args,kwargs)
       self.map_['source']= 'SDO'
       self.map_['provider'] ='LASP'
       self.map_['instrument'] = 'eve'
       self.map_['phyobs'] = 'irradiance'
    
