from sunpy.net.download  import Downloader
from sunpy.time import TimeRange
import datetime
from sunpy.net.vso.vso import Results
from sunpy.util import print_table
from sunpy.net.vso.attrs import Time

class queryrequestblock(object):
    """ Aims to provide user with additional information. Client.query returns
    container for queryrequestblock(s).It shall contain information about source,related URL.
    """
    def __init__(self, map_, url):
        self.source = map_.get('source', "Data not Available")
        self.provider = map_.get('provider', "Data not Available")
        self.phyobs = map_.get('phyobs', "Data not Available")
        self.instrument = map_.get('instrument', "Data not Available")
        self.url = url
        self.time = {}
        self.time['start'] = map_.get('Time_start', "Data not Available")
        self.time['end'] = map_.get('Time_end', "Data not available")

def iter_urls(map_, url_list):
    """Helper Function"""
    for aurl in url_list:
        tmp = queryrequestblock(map_, aurl)
        yield tmp


class queryresponse(list):
    """Returned by client.query.Attempt to ape QueryResponse object in vso module.
    """
    def __init__(self, lst):
    	
        super(queryresponse, self).__init__(lst)

    @classmethod
    def create(cls, map_, lst):
    	
	return cls(iter_urls(map_, lst))
    
    def time_range(self):
    	"""Returns the time-span query extends over"""
	return (datetime.date.strftime(
	        min(qrblock.time['start'] for qrblock in self), '%Y/%m/%d'),
		datetime.date.strftime(
		max(qrblock.time['end'] for qrblock in self), '%Y/%m/%d'))
    
    def show(self):
        """Presents data within container in a presentable manner"""
        table = [
	         [ 
		  qrblock.time['start'].strftime('%Y/%m/%d'),
		  qrblock.time['end'].strftime('%Y/%m/%d'),
		  qrblock.source,
		  qrblock.instrument,
		  qrblock.url
		 ] 
		  for qrblock in self
		]
	table.insert(0, ['Start time', 'End time', 'Source', 'Instrument', 'URL'])
	print print_table(table, colsep='  ', linesep='\n')



class GenericClient(object):    
    
    def __init__(self):
       self.map_ = {}

    def makeargs(self, *args, **kwargs):
       '''Convert Attribute in query to internal dictionary'''
       for elem in args:
           if issubclass(elem.__class__, Time):
               self.map_['TimeRange'] = TimeRange(elem.start, elem.end)
	       self.map_['Time_start'] = elem.start
	       self.map_['Time_end'] = elem.end
	   else:
	       try:
	           self.map_[elem.__class__.__name__] = elem.value
	       except Exception:
	           self.map_[elem.__class__.__name__] = None
       self._makeimap()

    
    def _get_url_for_timerange(cls, timerange, **kwargs):
       raise NotImplementedError

    def _get_url_for_date(cls, date, **kwargs):
       raise NotImplementedError
    
    @classmethod
    def _can_handle_query(cls, *query):
        raise NotImplementedError
    
    def query(self, *args, **kwargs):
         """
	 Input:
	 args: list of attributes

	 Output: queryresponse object.
	 """
	 GenericClient.makeargs(self, *args, **kwargs)
	 urls = self._get_url_for_timerange(self.map_.get('TimeRange'), **kwargs)
	 return queryresponse.create(self.map_, urls)

    
    def get(self, qres, **kwargs):
         """
	 Input:
	 qres : queryresponse object.

	 Output:
         vso.Results object.To wait for download to complete call .wait() on returned Results object.
	 """
	 urls = []
	 for qrblock in qres:
	     urls.append(qrblock.url)

         res = Results(lambda x: None, 0, lambda map_:self.link(map_))

	 dobj = Downloader(max_conn=len(urls), max_total=len(urls))
	 for aurl, ncall in list(zip(urls, map(lambda x:res.require([x]), urls))):
	     dobj.download(aurl, kwargs.get('Path',None), ncall, kwargs.get('ErrorBack', None))
         
	 return res

    def link(self, map_):
 	"""Helper Function"""
    	paths = []
	for k, v in map_.iteritems():
	    paths.append(map_[k]['path'])
	
	return paths
    
   
    def download_legacy(self, timerange, path=None, callback=None, errback=None):
	"""
	Aims to provide a simple way for downloading data by passing attribute usage.
	Input:
	timerange: Time-range over which to download data.
	path: Defaults to None in which case path used is one defined in sunpyrc file.
	callback: Function to be invoked at completion of download successfully.
	errorback: Function to be called when error is thrown during download.

	"""
        urls = EVEDownloader._get_url_for_timerange(timerange)
	dobj = Downloader(max_conn=len(urls), max_total=len(urls))
	for url in urls:
	    dobj.download(url, path, callback, errback)

