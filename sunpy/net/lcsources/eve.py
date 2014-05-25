import urlparse
from sunpy.net.vso.attrs import Time
from sunpy.net.download  import Downloader
from sunpy.time import TimeRange
import datetime
from sunpy.net.vso.vso import Results
from sunpy.util import print_table

class queryrequestblock(object):

    def __init__(self,url):
	self.source = 'SDO'
	self.provider = 'LASP'
	self.phyobs = 'Irradiance'
	self.instrument = 'EVE'
	self.url = url
	self.time = {}
	self.time['start'] = datetime.datetime.strptime(url.split('/')[-1].split('_')[0],'%Y%m%d')
	self.time['end'] = self.time['start'] + datetime.timedelta(days=1)


def iter_urls(url_list):

    for aurl in url_list:
          tmp = queryrequestblock(aurl)
	  print tmp
	  yield tmp


class queryresponse(list):

    def __init__(self,lst):
    	
	super(queryresponse,self).__init__(lst)

    @classmethod
    def create(cls,lst):
    	
        i= iter_urls(lst)
	return cls(iter_urls(lst))
    
    def num_records(self):
	
	return len(self)

    def time_range(self):
    	
	return (datetime.strftime(
	        min(qrblock.time['start'] for qrblock in self),'%Y/%m/%d'),
		datetime.strftime(
		max(qrblock.time['end'] for qrblock in self),'%Y/%m/%d'))
    
    def show(self):
    
        table = [
	         [ 
		  qrblock.time['start'].strftime('%Y/%m/%d'),
		  qrblock.time['end'].strftime('%Y/%m/%d'),
		  qrblock.source,
		  qrblock.instrument
		 ] 
		  for qrblock in self
		]
	table.insert(0,['Start time','End time','Source','Instrument'])
	print print_table(table, colsep='  ', linesep='\n')


class EVEDownloader(object):


    def __init__(self):
        
	self.map_={}

    @classmethod
    def _get_url_for_timerange(cls,timerange,**kwargs):
         days = timerange.get_days()
	 urls = []
	 for day in days:
	     urls.append(cls._get_url_for_date(day,**kwargs))
	 return urls

    @classmethod
    def _get_url_for_date(cls,date,**kwargs):
        
    #	if date < datetime.date(2010,1,1):
    #        raise ERROR (Define error class showing data not available for this date
        base_url = 'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/'
        return urlparse.urljoin(base_url, date.strftime('%Y/%Y%m%d') + '_EVE_L0CS_DIODES_1m.txt')
    
    def makeargs(self,*args,**kwargs):
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
         
	 self.makeargs(*args,**kwargs)
	 urls = EVEDownloader._get_url_for_timerange(self.map_.get('TimeRange'),**kwargs)
         return queryresponse.create(urls)

    
    def get(self,qres):
         
	 urls=[]
	 for qrblock in qres:
	     urls.append(qrblock.url)

         res = Results(lambda x: None,0,lambda map_:self.link(map_))

	 dobj = Downloader(max_conn=len(urls),max_total=len(urls))
	 for aurl,ncall in list(zip(urls,map(lambda x:res.require([x]),urls))):
	     dobj.download(aurl,self.map_.get('Path',None),ncall,self.map_.get('ErrorBack',None))
         
	 return res
         #add to database
	 #eve_add_to_database(conn_url)

    def link(self,map_):

    	paths = []
	for k,v in map_.iteritems():
	    paths.append(map_[k]['path'])
	
	return paths
    
   
    def download_legacy(self,timerange,path=None,callback=None,errback=None):
	
        urls = EVEDownloader._get_url_for_timerange(timerange)
	print urls
	print len(urls)
	dobj = Downloader(max_conn=len(urls),max_total=len(urls))
	for url in urls:
	    dobj.download(url,path,callback,errback)
     

     



	 


