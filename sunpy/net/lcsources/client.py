from sunpy.net.download  import Downloader
from sunpy.time import TimeRange
import datetime
from sunpy.net.vso.vso import Results
from sunpy.util import print_table
from sunpy.net.vso.attrs import Time

class queryrequestblock(object):

    def __init__(self,map_,url):
	self.source = map_.get('source',"Data not Available")
	self.provider = map_.get('provider',"Data not Available")
	self.phyobs = map_.get('phyobs',"Data not Available")
	self.instrument = map_.get('instrument',"Data not Available")
	self.url = url
	self.time = {}
	self.time['start'] = datetime.datetime.strptime(url.split('/')[-1].split('_')[0],'%Y%m%d')
	self.time['end'] = self.time['start'] + datetime.timedelta(days=1)

def iter_urls(map_,url_list):

    for aurl in url_list:
          tmp = queryrequestblock(map_,aurl)
	  yield tmp


class queryresponse(list):

    def __init__(self,lst):
    	
	super(queryresponse,self).__init__(lst)

    @classmethod
    def create(cls,map_,lst):
    	
	return cls(iter_urls(map_,lst))
    
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



class GenericClient(object):    
    
    def __init__(self):
       self.map_ = {}

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
    
    def _get_url_for_timerange(cls,timerange,**kwargs):
       
       raise NotImplementedError

    def _get_url_for_date(cls,date,**kwargs):
   	
       raise NotImplementedError

    
    def query(self,*args,**kwargs):
         
	 GenericClient.makeargs(self,*args,**kwargs)
	 urls = self._get_url_for_timerange(self.map_.get('TimeRange'),**kwargs)
	 return queryresponse.create(self.map_,urls)

    
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
	dobj = Downloader(max_conn=len(urls),max_total=len(urls))
	for url in urls:
	    dobj.download(url,path,callback,errback)

