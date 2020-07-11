from pathlib import Path
from collections import OrderedDict

import astropy.table
import astropy.units as u
from astropy.time import TimeDelta

import sunpy
from sunpy import config
from sunpy.net._attrs import Time, Wavelength
from sunpy.net.base_client import BaseClient, BaseQueryResponse
from sunpy.time import TimeRange, parse_time
from sunpy.util.parfive_helpers import Downloader
from sunpy.util.scraper import Scraper

TIME_FORMAT = config.get("general", "time_format")

__all__ = ['QueryResponse', 'GenericClient']


class QueryResponseBlock(OrderedDict):
    """
    Represents each row for the QueryResponse table.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for k in self:
            self.__dict__[k.lower()] = self[k]


class QueryResponse(BaseQueryResponse):
    """
    A container for files metadata returned by
    searches from Dataretriver clients.
    """

    def __init__(self, lst, client=None):
        super().__init__()
        datalist = []
        for datablock in lst:
            datalist.append(QueryResponseBlock(datablock))
        self._data = datalist
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

    def time_range(self):
        """
        Returns the time-span for which records are available.
        """
        return TimeRange(min(qrblock['Time'].start for qrblock in self),
                         max(qrblock['Time'].end for qrblock in self))

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

        # finding column names to be shown in the response table.
        colnames = []
        for colname in self._data[0].keys():
            if colname != 'url' and colname != 'Time':
                colnames.append(colname)
        columns = OrderedDict(((col, [])) for col in colnames)

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
    baseurl = None
    pattern = None

    @classmethod
    def _get_match_dict(cls, *args, **kwargs):
        """
        Returns a dictionary to validate the metadata of searched files using
        query and registered values for Attrs for the respective client.

        Parameters
        ----------
        \\*args: `tuple`
            `sunpy.net.attrs` objects representing the query.
        \\*\\*kwargs: `dict`
             Any extra keywords to refine the search.

        Returns
        -------
        matchdict: `dict`
            A dictionary having a `list` of all possible Attr values
            corresponding to an Attr.
        """
        a = cls.register_values()
        matchdict = {}
        for i in a.keys():
            attrname = i.__name__
            matchdict[attrname] = []
            for val, desc in a[i]:
                matchdict[attrname].append(val)
        for elem in args:
            if isinstance(elem, Time):
                timerange = TimeRange(elem.start, elem.end)
                matchdict['Time'] = timerange
            elif hasattr(elem, 'value'):
                matchdict[elem.__class__.__name__] = [str(elem.value).lower()]
            elif isinstance(elem, Wavelength):
                matchdict['Wavelength'] = elem
            else:
                raise ValueError("GenericClient can not add {} to the map_ dictionary to pass to the Client.".format(elem.__class__.__name__))
        for k in kwargs:
            matchdict[k] = [kwargs[k]]
            if isinstance(kwargs[k], str):
                matchdict[k] = [kwargs[k].lower()]
        return matchdict

    @classmethod
    def pre_search_hook(cls, *args, **kwargs):
        """
        Helper function to return the baseurl, pattern and matchdict
        for the client required by `search()` before using the scraper.
        """
        matchdict = cls._get_match_dict(*args, **kwargs)
        return cls.baseurl, cls.pattern, matchdict

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
        for key in adict:
            all_vals = [i[0].lower() for i in adict[key]]
            for x in query:
                if isinstance(x, key) and str(x.value).lower() not in all_vals:
                    return False
        return True

    def post_search_hook(self, exdict, matchdict):
        """
        Helper function used after `search()` which makes the
        extracted metadata representable in a query response table.

        Parameters
        ----------
        exdict: `dict`
            Represents metadata extracted from files.
        matchdict: `dict`
            Contains attr values accessed from `register_values()`
            and the search query itself.

        Returns
        -------
        map_: `~collections.OrderedDict`
            An Ordered Dictionary which is used by `QueryResponse`
            to show results.
        """
        map_ = OrderedDict()
        almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
        if 'month' in exdict and 'day' in exdict:
            start = parse_time("{}/{}/{}".format(exdict['year'], exdict['month'], exdict['day']))
            end = start + almost_day
        elif 'year' in exdict:
            start = parse_time("{}/{}/{}".format(exdict['year'], 1, 1))
            end = parse_time("{}/{}/{}".format(exdict['year'], 12, 31)) + almost_day
        map_['Time'] = TimeRange(start, end)
        map_['Start Time'] = start.strftime(TIME_FORMAT)
        map_['End Time'] = end.strftime(TIME_FORMAT)
        map_['Instrument'] = matchdict['Instrument'][0].upper()
        if 'Physobs' in matchdict:
            map_['Physobs'] = matchdict['Physobs'][0]
        map_['Source'] = matchdict['Source'][0].upper()
        map_['Provider'] = matchdict['Provider'][0].upper()
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

            temp_dict = qres.blocks[i].copy()
            temp_dict['file'] = str(filename)
            fname = fname.expanduser()
            fname = Path(str(fname).format(**temp_dict))

            paths.append(fname)

        return paths

    def search(self, *args, **kwargs):
        """
        Query this client for a list of results.

        Parameters
        ----------
        \\*args: `tuple`
            `sunpy.net.attrs` objects representing the query.
        \\*\\*kwargs: `dict`
             Any extra keywords to refine the search.

        Returns
        -------
        A `QueryResponse` instance containing the query result.
        """
        baseurl, pattern, matchdict = self.pre_search_hook(*args, **kwargs)
        scraper = Scraper(baseurl, regex=True)
        filesmeta = scraper._extract_files_meta(matchdict['Time'], extractor=pattern, matcher=matchdict)
        metalist = []
        for i in filesmeta:
            map_ = self.post_search_hook(i, matchdict)
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

        urls = [qrblock['url'] for qrblock in qres]

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
