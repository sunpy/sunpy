#Author :Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#This module was developed under funding provided by
#Google Summer of Code 2014

import copy
import os
import datetime
from collections import OrderedDict
from functools import partial

import astropy.table

import sunpy
from sunpy.extern import six
from sunpy.time import TimeRange
from sunpy.util import replacement_filename
from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

from ..download import Downloader, Results
from ..vso.attrs import Time

__all__ = ['QueryResponse', 'GenericClient']



def simple_path(path, sock, url):
    return path

class QueryResponseBlock(object):
    """
    Represents url, source along with other information
    """

    def __init__(self, map0, url):
        """
        Parameters
        ----------
        map0 : Dict with relevant information
        url  : Uniform Resource Locator
        """
        self.source = map0.get('source', "Data not Available")
        self.provider = map0.get('provider', "Data not Available")
        self.phyobs = map0.get('phyobs', "Data not Available")
        self.instrument = map0.get('instrument', "Data not Available")
        self.url = url
        self.time = TimeRange(map0.get('Time_start'), map0.get('Time_end'))


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

    def __repr__(self):
        return repr(self._build_table())

    def __str__(self):
        return str(self._build_table())

    def _repr_html_(self):
        return self._build_table()._repr_html_()

    def _build_table(self):
        columns = OrderedDict((('Start Time', []), ('End Time', []), (
            'Source', []), ('Instrument', [])))
        for i, qrblock in enumerate(self):
            columns['Start Time'].append((qrblock.time.start.date(
            ) + datetime.timedelta(days=i)).strftime(TIME_FORMAT))
            columns['End Time'].append((qrblock.time.start.date(
            ) + datetime.timedelta(days=i+1)).strftime(TIME_FORMAT))
            columns['Source'].append(qrblock.source)
            columns['Instrument'].append(qrblock.instrument)

        return astropy.table.Table(columns)


class GenericClient(object):
    """
    Base class for simple web clients for the data retriever module. This class
    is mainly designed for downloading data from FTP and HTTP type data
    sources, although should in theory be general enough to get data from any
    web service.

    This class has two user facing methods
    `~sunpy.net.dataretriever.client.GenericClient.query` and
    `~sunpy.net.dataretriever.client.GenericClient.get` the former generates a
    set of results for files available through the service the client is
    querying and the latter downloads that data.

    The `~sunpy.net.dataretriever.client.GenericClient.query` method takes a
    set of `sunpy.net.attrs` objects and then converts these into a call to
    `~sunpy.net.dataretriever.client.GenericClient._get_url_for_timerange`. It
    does this through the `map\_` dictionary which represents the
    `~sunpy.net.attrs` objects as a dictionary.
    """

    def __init__(self):
        self.map_ = {}

    def _makeargs(self, *args, **kwargs):
        """
        Construct the `map\_` internal representation of the query.

        This `map\_` dictionary is passed through to the
        `_get_url_for_timerange` method to get the URL results.

        Parameters
        ----------
        \*args: `tuple`
            The query attributes.

        \*\*kwargs: `dict`
            None.
        """
        for elem in args:
            if issubclass(elem.__class__, Time):
                self.map_['TimeRange'] = TimeRange(elem.start, elem.end)
                self.map_['Time_start'] = elem.start
                self.map_['Time_end'] = elem.end
            else:
                try:
                    self.map_[elem.__class__.__name__.lower()] = elem.value
                except Exception:
                    self.map_[elem.__class__.__name__.lower()] = None
        self._makeimap()

    def _get_url_for_timerange(cls, timerange, **kwargs):
        """
        Method which generates URL results from a timerange and the `map\_`
        dictionary.

        Parameters
        ----------
        timerange: `sunpy.time.TimeRange`
             The timerange to extract the URLs for.

        \*\*kwargs: `dict`
             Any extra keywords to refine the search. Generated from the
             attributes passed to
             `~sunpy.net.dataretriever.client.GenericClient.query`.
        """
        raise NotImplementedError

    def _makeimap(self):
        """
        Add client specific information to the _map dict.

        Normally this is extra metadata which is not downloaded, but known
        a priori.
        """
        raise NotImplementedError

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Method the
        `sunpy.net.dataretriever.downloader_factory.UnifiedDownloaderFactory`
        class uses to dispatch queries to this Client.
        """
        raise NotImplementedError

    def query(self, *args, **kwargs):
        """
        Query this client for a list of results.

        Parameters
        ----------
        \*args: `tuple`
            `sunpy.net.attrs` objects representing the query.
        """
        GenericClient._makeargs(self, *args, **kwargs)

        kwergs = copy.copy(self.map_)
        kwergs.update(kwargs)
        urls = self._get_url_for_timerange(
            self.map_.get('TimeRange'), **kwergs)
        return QueryResponse.create(self.map_, urls)

    def get(self, qres, path=None, error_callback=None, **kwargs):
        """
        Download a set of results.

        Parameters
        ----------
        qres : `~sunpy.net.dataretriever.QueryResponse`
            Results to download.

        Returns
        -------
        Results Object
        """
        details = qres.client.map_

        urls = []
        for qrblock in qres:
            urls.append(qrblock.url)

        filenames = []
        for url in urls:
            filenames.append(url.split('/')[-1])

        # Create function to compute the filepath to download to if not set
        default_dir = sunpy.config.get("downloads", "download_dir")

        paths = []
        for filename in filenames:
            if path is None:
                path = os.path.join(default_dir, '{file}')
            elif isinstance(path, six.string_types) and '{file}' not in path:
                path = os.path.join(path, '{file}')

            temp_dict = details.copy()
            temp_dict['file'] = filename
            path  = path.format(**temp_dict)
            path = os.path.expanduser(path)

            if os.path.exists(path):
                path = replacement_filename(path)

            path = partial(simple_path, path)

            paths.append(path)

        res = Results(lambda x: None, 0, lambda map_: self._link(map_))

        dobj = Downloader(max_conn=len(urls), max_total=len(urls))

        # We cast to list here in list(zip... to force execution of 
        # res.require([x]) at the start of the loop.
        for aurl, ncall, path in list(zip(urls, map(lambda x: res.require([x]),
                                              urls), paths)):
            dobj.download(aurl, path, ncall, error_callback)

        return res

    def _link(self, map_):
        """Helper Function"""
        paths = []
        for k, v in map_.iteritems():
            paths.append(map_[k]['path'])
        return paths
