from pathlib import Path
from collections import OrderedDict

import astropy.table
import astropy.units as u
from astropy.time import TimeDelta

import sunpy
from sunpy import config
from sunpy.net._attrs import Time
from sunpy.net.base_client import BaseClient, BaseQueryResponse
from sunpy.time import TimeRange, parse_time
from sunpy.util.parfive_helpers import Downloader
from sunpy.util.scraper import Scraper

TIME_FORMAT = config.get("general", "time_format")

__all__ = ['QueryResponse', 'GenericClient']


class QueryResponse(BaseQueryResponse):

    def __init__(self, lst, client=None):
        super().__init__()
        self._data = lst
        self._client = client

    @property
    def blocks(self):
        return self._data

    @property
    def client(self):
        return self._client

    @client.setter
    def client(self, client):
        self._client = client

    def __len__(self):
        return len(self._data)

    def __getitem__(self, item):
        # Always index so a list comes back
        if isinstance(item, int):
            item = slice(item, item+1)
        return type(self)(self._data[item], client=self.client)

    def __iter__(self):
        for block in self._data:
            yield block

    def response_block_properties(self):
        """
        Returns a set of class attributes on all the response blocks.
        """
        s = {a if not a.startswith('_') else None for a in dir(self[0])}
        for resp in self[1:]:
            s = s.intersection({a if not a.startswith('_') else None for a in dir(resp)})

        s.remove(None)
        return s

    def build_table(self):
        if len(self._data) == 0:
            return astropy.table.Table()
        meta0 = self._data[0]
        meta0.pop('url', None)
        columns = OrderedDict(((col, [])) for col in meta0.keys())
        for qrblock in self:
            for colname in columns.keys():
                columns[colname].append(qrblock[colname])
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
    """

    @classmethod
    def pre_hook(cls, *args, **kwargs):
        """
        """
        a = cls.register_values()
        d = {}
        for i in a.keys():
            attrname = i.__name__
            d[attrname] = []
            for val, desc in a[i]:
                d[attrname].append(val)
        for elem in args:
            if isinstance(elem, Time):
                timerange = TimeRange(elem.start, elem.end)
                d['timerange'] = timerange
            elif hasattr(elem, 'value'):
                d[elem.__class__.__name__] = [str(elem.value)]
            else:
                raise ValueError("GenericClient can not add {} to the map_ dictionary to pass to the Client.".format(elem.__class__.__name__))
        for k in kwargs:
            d[k] = [kwargs[k]]
        return cls.baseurl, cls.pattern, d

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Method the
        `sunpy.net.fido_factory.UnifiedDownloaderFactory`
        class uses to dispatch queries to this Client.
        """
        from sunpy.net import attrs as a

        required = {a.Time, a.Instrument}
        adict = cls.register_values()
        optional = {i for i in adict.keys()} - required
        if hasattr(cls, 'required'):
            required = cls.required
        if hasattr(cls, 'optional'):
            optional = cls.optional
        if not cls.check_attr_types_in_query(query, required, optional):
            return False
        for x in query:
            if isinstance(x, a.Instrument) and x.value.lower() == adict[a.Instrument][0][0].lower():
                return True
        return False

    def post_hook(self, exdict, matchdict):
        """
        """
        map_ = OrderedDict()
        start = parse_time("{}/{}/{}".format(exdict['year'], exdict['month'], exdict['day']))
        end = start + TimeDelta(1 * u.day - 1 * u.millisecond)
        map_['Start Time'] = start.strftime(TIME_FORMAT)
        map_['End Time'] = end.strftime(TIME_FORMAT)
        map_['Instrument'] = matchdict['Instrument'][0]
        if 'Physobs' in matchdict:
            map_['Phsyobs'] = matchdict['Physobs'][0]
        map_['Source'] = matchdict['Source'][0]
        map_['Provider'] = matchdict['Provider'][0]
        for k in exdict:
            if k not in ['year', 'month', 'day']:
                map_[k] = exdict[k]
        return map_

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

            temp_dict = qres.blocks[i]._map.copy()
            temp_dict['file'] = str(filename)
            fname = fname.expanduser()
            fname = Path(str(fname).format(**temp_dict))

            paths.append(fname)

        return paths

    def search(self, *args, **kwargs):
        """
        Query this client for a list of results.
        """
        baseurl, pattern, matchdict = self.pre_hook(*args, **kwargs)
        scraper = Scraper(baseurl, regex=True)
        filesmeta = scraper._extract_files_meta(matchdict['timerange'], extractor=pattern, matcher=matchdict)
        metalist = []
        for i in filesmeta:
            map_ = self.post_hook(i, matchdict)
            metalist.append(map_)
        return QueryResponse(metalist, client=self)

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
        overwrite : `bool` or `str`, optional
            Determine how to handle downloading if a file already exists with the
            same name. If `False` the file download will be skipped and the path
            returned to the existing file, if `True` the file will be downloaded
            and the existing file will be overwritten, if `'unique'` the filename
            will be modified to be unique.
        progress : `bool`, optional
            If `True` show a progress bar showing how many of the total files
            have been downloaded. If `False`, no progress bar will be shown.
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

        urls = [qrblock.url for qrblock in qres.blocks]

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
