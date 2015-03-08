#Author :Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#This module was developed under funding provided by
#Google Summer of Code 2014

import datetime

from sunpy.time import TimeRange
from sunpy.util import print_table

from sunpy.net.download import Downloader
from sunpy.net.vso.vso import Results
from sunpy.net.vso.attrs import Time

__all__ = ['QueryResponse', 'GenericClient']


class QueryResponseBlock(object):
    """
    Represents url, source along with other information
    """
    def __init__(self, map_, url):
        """
        Parameters
        ----------
        map_ : Dict with relevant information
        url  : Uniform Resource Locator
        """
        self.source = map_.get('source', "Data not Available")
        self.provider = map_.get('provider', "Data not Available")
        self.phyobs = map_.get('phyobs', "Data not Available")
        self.instrument = map_.get('instrument', "Data not Available")
        self.url = url
        self.time = TimeRange(map_.get('Time_start'), map_.get('Time_end'))


def iter_urls(map_, url_list):
    """Helper Function"""
    for aurl in url_list:
        tmp = QueryResponseBlock(map_, aurl)
        yield tmp


class QueryResponse(list):
    """
    Container of QueryResponseBlocks
    """
    def __init__(self, lst):

        super(QueryResponse, self).__init__(lst)

    @classmethod
    def create(cls, map_, lst):

        return cls(iter_urls(map_, lst))

    def time_range(self):
        """
        Returns the time-span for which records are available
        """
        return (datetime.date.strftime(
                min(qrblock.time.start for qrblock in self), '%Y/%m/%d'),
                datetime.date.strftime(
                max(qrblock.time.end for qrblock in self), '%Y/%m/%d'))

    def __str__(self):
        """Presents data within container in a presentable manner"""

        table = [
                 [
                  (qrblock.time.t1.date() + datetime.timedelta(days=i)).strftime('%Y/%m/%d'),
                  (qrblock.time.t2.date() + datetime.timedelta(days=i)).strftime('%Y/%m/%d'),
                  qrblock.source,
                  qrblock.instrument,
                  qrblock.url
                 ]
                  for i,qrblock in enumerate(self)
                ]
        table.insert(0, ['----------', '--------', '------', '----------', '---'])
        table.insert(0, ['Start time', 'End time', 'Source', 'Instrument', 'URL'])
        return print_table(table, colsep='  ', linesep='\n')


class GenericClient(object):

    def __init__(self):
        self.map_ = {}

    def _makeargs(self, *args, **kwargs):
        '''Map attributes in the query to internal dictionary'''
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
        Query the web service of the source for urls pertaining to incoming arguements.
        """
        GenericClient._makeargs(self, *args, **kwargs)
        urls = self._get_url_for_timerange(self.map_.get('TimeRange'), **kwargs)
        return QueryResponse.create(self.map_, urls)


    def get(self, qres, **kwargs):
        """
        Parameters
        ----------
        qres : QueryResponse object

        Returns
        -------
        Results Object
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
        Download required data using keyword arguments.

        Parameters
        ----------
        timerange: Time-range over which to download data.
        path: Defaults to None in which case path used is one defined in sunpyrc file.
        callback: Function to be invoked at completion of download successfully.
        errorback: Function to be called when error is thrown during download.

        Examples
        --------
        >>> import sunpy.net.unifieddownloader.sources.eve as eve
        >>> cl = eve.EVEClient()
        >>> cl.download_legacy(Time('2012/2/2','2012/2/3'))
        """
        urls = self._get_url_for_timerange(timerange)
        dobj = Downloader(max_conn=len(urls), max_total=len(urls))
        for url in urls:
            dobj.download(url, path, callback, errback)
