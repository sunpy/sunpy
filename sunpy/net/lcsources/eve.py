import urlparse
from sunpy.net.vso.attrs import Time
from sunpy.net.download  import Downloader



class EVEDownloader(object):

    def _get_url_from_timerange(self,timerange,**kwargs):
         days = timerange.get_days()
	 urls = []
	 for day in days:
	     urls.append(cls._get_url_for_date(day,**kwargs))
	 return urls

    def _get_url_for_date(self,date,**kwargs):
        
    #	if date < datetime.date(2010,1,1):
    #        raise ERROR (Define error class showing data not available for this date
        base_url = 'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/'
        return urlparse.urljoin(base_url, date.strftime('%Y/%Y%m%d') + '_EVE_L0CS_DIODES_1m.txt')
    
    @staticmethod
    def makeargs(map_,*args,**kwargs):
       '''map_:Dict'''
       for elem in args:
           if issubclass(elem.__class__,Time):
               map_['TimeRange']=TimeRange(elem.start,elem.end)
	   else:
	       try:
	           map_[attr.__class__.__name__]=attr.value
	       except Exception:
	           map_[attr.__class__.__name__]=None
     
    def download(self,*args,**kwargs):
         
         map_ = {}
	 makeargs(map_,*args,**kwargs)
	 urls = EVEDownloader._get_url_for_timerange(map_.get('TimeRange'),**kwargs)
	 dobj = Downloader(max_conn=len(urls),max_total=len(urls))
	 for url in urls:
	     dobj.download(url,map_.get('Path',None),map_.get('CallBack',None),map_.get('ErrorBack',None))
         
         #add to database
	 #eve_add_to_database(conn_url)

     
    @staticmethod
    def eve_add_to_database(conn_url):
	  pass
    
   
    def download_legacy(self,timerange,path=None,callback=None,errback=None):

	 print dir(EVEDownloader)
         urls = self._get_url_for_timerange(timerange)
	 dobj = Downloader(max_conn=len(urls),max_total=len(urls))
	 for url, in urls:
	     dobj.download(url,path,callback,errback)
     

     



	 


