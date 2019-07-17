# -*- coding: utf-8 -*-
import copy
from pathlib import Path
from collections import OrderedDict, namedtuple

import numpy as np
from parfive import Downloader

import astropy.table
import astropy.units as u

import parfive

import sunpy
from sunpy import config
from sunpy.net.base_client import BaseClient
from sunpy.net.vso.attrs import Time, Wavelength, _Range
from sunpy.time import TimeRange

TIME_FORMAT = config.get("general", "time_format")

__all__ = ['QueryResponse', 'GenericClient']


class QueryResponseBlock(object):
    """
    Represents url, source along with other information
    """

    def __init__(self, map0, url, time=None):
        """
        Parameters
        ----------
        map0 : Dict with relevant information
        url  : Uniform Resource Locator
        """
        self._map = map0
        self.source = map0.get('source', "Data not Available")
        self.provider = map0.get('provider', "Data not Available")
        self.physobs = map0.get('physobs', "Data not Available")
        self.instrument = map0.get('instrument', "Data not Available")
        self.url = url
        self.time = TimeRange(map0.get('Time_start'),
                              map0.get('Time_end')) if time is None else time
        self.wave = map0.get('wavelength', np.NaN)


def iter_urls(amap, url_list, time):
    """Helper Function"""
    for aurl, t in zip(url_list, time):
        tmp = QueryResponseBlock(amap, aurl, t)
        yield tmp


class QueryResponse(list):
    """
    Container of QueryResponseBlocks
    """

    def __init__(self, lst):
        super(QueryResponse, self).__init__(lst)

    @classmethod
    def create(cls, amap, lst, time=None):
        if time is None:
            time = [None] * len(lst)
        return cls(iter_urls(amap, lst, time))

    def time_range(self):
        """
        Returns the time-span for which records are available
        """
        return TimeRange(min(qrblock.time.start for qrblock in self),
                         max(qrblock.time.end for qrblock in self))

    def response_block_properties(self):
        """
        Returns a set of class attributes on all the response blocks.
        """
        s = {a if not a.startswith('_') else None for a in dir(self[0])}
        for resp in self[1:]:
            s = s.intersection({a if not a.startswith('_') else None for a in dir(resp)})

        s.remove(None)
        return s

    def __repr__(self):
        return repr(type(self)) + repr(self._build_table())

    def __str__(self):
        return str(self._build_table())

    def _repr_html_(self):
        return self._build_table()._repr_html_()

    def _build_table(self):
        columns = OrderedDict((('Start Time', []), ('End Time', []),
                               ('Source', []), ('Instrument', []),
                               ('Wavelength', [])))
        for i, qrblock in enumerate(self):
            columns['Start Time'].append(
                (qrblock.time.start).strftime(TIME_FORMAT))
            columns['End Time'].append(
                (qrblock.time.end).strftime(TIME_FORMAT))
            columns['Source'].append(qrblock.source)
            columns['Instrument'].append(qrblock.instrument)
            columns['Wavelength'].append(str(u.Quantity(qrblock.wave)))

        return astropy.table.Table(columns)


class GenericClient(BaseClient):
    """
    Base class for simple web clients for the data retriever module. This class
    is mainly designed for downloading data from FTP and HTTP type data
    sources, although should in theory be general enough to get data from any
    web service.

    This class has two user facing methods
    `~sunpy.net.dataretriever.client.GenericClient.search` and
    `~sunpy.net.dataretriever.client.GenericClient.fetch` the former generates a
    set of results for files available through the service the client is
    querying and the latter downloads that data.

    The `~sunpy.net.dataretriever.client.GenericClient.search` method takes a
    set of `sunpy.net.attrs` objects and then converts these into a call to
    `~sunpy.net.dataretriever.client.GenericClient._get_url_for_timerange`. It
    does this through the `map\\_` dictionary which represents the
    `~sunpy.net.attrs` objects as a dictionary.
    """

    def __init__(self):
        self.map_ = {}

    def _makeargs(self, *args):
        """
        Construct the `map\\_` internal representation of the query.

        This `map\\_` dictionary is passed through to the
        `_get_url_for_timerange` method to get the URL results.

        Parameters
        ----------
        \\*args: `tuple`
            The query attributes.

        """
        for elem in args:
            if isinstance(elem, Time):
                self.map_['TimeRange'] = TimeRange(elem.start, elem.end)
                self.map_['Time_start'] = elem.start
                self.map_['Time_end'] = elem.end
            elif isinstance(elem, _Range):
                a_min = elem.min
                a_max = elem.max
                if a_min == a_max:
                    self.map_[elem.__class__.__name__.lower()] = a_min
                else:
                    if isinstance(elem, Wavelength):
                        prefix = 'wave'
                    else:
                        prefix = ''
                    minmax = namedtuple("minmax", "{0}min {0}max".format(prefix))
                    self.map_[elem.__class__.__name__.lower()] = minmax(a_min, a_max)
            else:
                if hasattr(elem, 'value'):
                    self.map_[elem.__class__.__name__.lower()] = elem.value
                else:
                    # This will only get hit if the attr is something like
                    # Extent, which is a unique subclass of Attr. Currently no
                    # unidown Clients support this, so we skip this line.
                    # Anything that hits this will require special code to
                    # convert it into the map_ dict.
                    raise ValueError(
                        "GenericClient can not add {} to the map_ dictionary to pass "
                        "to the Client.".format(elem.__class__.__name__))  # pragma: no cover
        self._makeimap()

    @classmethod
    def _get_url_for_timerange(cls, timerange, **kwargs):
        """
        Method which generates URL results from a timerange and the `map\\_`
        dictionary.

        Parameters
        ----------
        timerange: `sunpy.time.TimeRange`
             The timerange to extract the URLs for.

        \\*\\*kwargs: `dict`
             Any extra keywords to refine the search. Generated from the
             attributes passed to
             `~sunpy.net.dataretriever.client.GenericClient.search`.
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
        `sunpy.net.fido_factory.UnifiedDownloaderFactory`
        class uses to dispatch queries to this Client.
        """
        raise NotImplementedError

    def _get_full_filenames(self, qres, filenames, path):
        """
        Download a set of results.

        Parameters
        ----------
        qres : `~sunpy.net.dataretriever.QueryResponse`
            Results to download.

        filenames : list
            List of base filenames (ex - "xyz.txt")

        path : str
            Path to download files to

        Returns
        -------
        List of full pathnames for each file (download_directory + filename)
        """
        # Create function to compute the filepath to download to if not set
        default_dir = Path(sunpy.config.get("downloads", "download_dir"))

        paths = []
        for i, filename in enumerate(filenames):
            fname = Path(filename)
            if path is None:
                fname = default_dir / '{file}'
            elif '{file}' not in str(path):
                fname = path / '{file}'

            temp_dict = qres[i]._map.copy()
            temp_dict['file'] = str(filename)
            fname = fname.expanduser()
            fname = Path(str(fname).format(**temp_dict))

            paths.append(fname)

        return paths

    def _get_time_for_url(self, urls):
        """
        This method allows clients to customise the timerange displayed for
        each URL.

        It should return a sunpy.time.TimeRange object per URL.
        """
        return NotImplemented

    def search(self, *args, **kwargs):
        """
        Query this client for a list of results.

        Parameters
        ----------
        \\*args: `tuple`
            `sunpy.net.attrs` objects representing the query.
        """
        GenericClient._makeargs(self, *args, **kwargs)

        kwergs = copy.copy(self.map_)
        kwergs.update(kwargs)
        urls = self._get_url_for_timerange(
            self.map_.get('TimeRange'), **kwergs)
        if urls:
            times = self._get_time_for_url(urls)
            if times and times is not NotImplemented:
                return QueryResponse.create(self.map_, urls, times)
        return QueryResponse.create(self.map_, urls)

    def fetch(self, qres, path=None, overwrite=False,
              progress=True, downloader=None, wait=True):
        """
        Download a set of results.

        Parameters
        ----------
        qres : `~sunpy.net.dataretriever.QueryResponse`
            Results to download.

        path : `str` or `pathlib.Path`, optional
            Path to the download directory, or file template including the
            ``{file}`` string which will be replaced with the filename.

        progress : `bool`, optional
            If `True` show a progress bar showing how many of the total files
            have been downloaded. If `False`, no progress bar will be shown.

        overwrite : `bool` or `str`, optional
            Determine how to handle downloading if a file already exists with the
            same name. If `False` the file download will be skipped and the path
            returned to the existing file, if `True` the file will be downloaded
            and the existing file will be overwritten, if `'unique'` the filename
            will be modified to be unique.

        downloader : `parfive.Downloader`, optional
            The download manager to use.

        wait : `bool`, optional
           If `False` ``downloader.download()`` will not be called. Only has
           any effect if `downloader` is not `None`.

        Returns
        -------

        results: `parfive.Results`

        """
        if path is not None:
            path = Path(path)

        urls = [qrblock.url for qrblock in qres]

        filenames = [url.split('/')[-1] for url in urls]

        paths = self._get_full_filenames(qres, filenames, path)

        dl_set = True
        if not downloader:
            dl_set = False
            downloader = Downloader(progress=progress, overwrite=overwrite)

        for url, filename in zip(urls, paths):
            downloader.enqueue_file(url, filename=filename)

        if dl_set and not wait:
            return

        return downloader.download()

    def _link(self, map_):
        """Helper Function"""
        paths = []
        for k, v in map_.items():
            paths.append(map_[k]['path'])
        return paths
