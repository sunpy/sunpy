import urlparse
from sunpy.net.vso.attrs import Time
from sunpy.net.download  import Downloader
from sunpy.time import TimeRange
import datetime



class queryrequestblock(object):

    def __init(self,url):
	self.source = 'SDO'
	self.provider = 'LASP'
	self.phyobs = 'Irradiance'
	self.instrument = 'EVE'
	self.url = url
	self.time = {}
	self.time['start'] = datetime(url.split('/')[-1].split('_')[0]).strpftime('%Y%m%d')
	self.time['end '] = self.time['start'] + datetime.timedelta(days=1)


def iter_urls(url_list):

    for aurl in url_list:
          tmp = queryrequestblock(url)
	  yield tmp







class queryresponse(list):

    def __init__(self,lst):
    	
	super(queryresponse,self).__init__(lst)

    @classmethod
    def create(cls,lst):
    	
	return cls(iter_urls(lst))
    
    def num_records(self):
	
	return self

    def time_range(self):
    	
	return (datetime.strftime(
	        min(qrblock.time['start'] for qrblock in self),'%Y/%m/%d'),
		datetime.strftime(
		max(qrblock.time['end'] for qrblock in self),'%Y/%m/%d'))
    
    def show():
    
        table = [
	          qrblock.time['start'].strftime('%Y/%m/%d'),
		  qrblock.time['end'].strftime('%Y/%m/%d'),
		  qrblock.source,
		  qrblock.instrument,
		]
	table.insert(0,[])
	table.insert(0,['Start time','End time','Source','Instrument'])
	print print_table(table, colsep='  ', linesep-'\n')




#TODO: add wait() method which blocks the download


class EVEDownloader(object):


    def __init__(self):
        
	self.map_={}

    @classmethod
    def _get_url_for_timerange(cls,timerange,**kwargs):
         days = timerange.get_days()
	 urls = []
	 print days
	 for day in days:
	     urls.append(cls._get_url_for_date(day,**kwargs))
	 print urls
	 return urls

    @classmethod
    def _get_url_for_date(cls,date,**kwargs):
        
    #	if date < datetime.date(2010,1,1):
    #        raise ERROR (Define error class showing data not available for this date
        base_url = 'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/'
        return urlparse.urljoin(base_url, date.strftime('%Y/%Y%m%d') + '_EVE_L0CS_DIODES_1m.txt')
    
    @staticmethod
    def makeargs(*args,**kwargs):
       '''map_:Dict'''
       for elem in args:
           if issubclass(elem.__class__,Time):
               self.map_['TimeRange']=TimeRange(elem.start,elem.end)
	   else:
	       try:
	           self.map_[elem.__class__.__name__]=elem.value
	       except Exception:
	           self.map_[elem.__class__.__name__]=None
     
    def query(self,*args,**kwargs):
         
	 EVEDownloader.makeargs(*args,**kwargs)
	 urls = EVEDownloader._get_url_for_timerange(self.map_.get('TimeRange'),**kwargs)
         return queryresponse.create(urls)


    def get(self,qres):
         
	 urls=[]
	 for qrblock in qres:
	     urls.append(qrblock.url)

	 dobj = Downloader(max_conn=len(urls),max_total=len(urls))
	 for url in urls:
	     dobj.download(url,self.map_.get('Path',None),self.map_.get('CallBack',None),self.map_.get('ErrorBack',None))
         
         #add to database
	 #eve_add_to_database(conn_url)

     
    
   
    def download_legacy(self,timerange,path=None,callback=None,errback=None):
	
	print "hello"
        urls = EVEDownloader._get_url_for_timerange(timerange)
	print urls
	print len(urls)
	dobj = Downloader(max_conn=len(urls),max_total=len(urls))
	for url in urls:
	    dobj.download(url,path,callback,errback)
     

     



	 


