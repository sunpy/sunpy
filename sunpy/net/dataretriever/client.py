from pathlib import Path
from collections import OrderedDict

import numpy as np

import sunpy
from sunpy import config
from sunpy.net import attrs as a
from sunpy.net.attr import SimpleAttr
from sunpy.net.base_client import BaseClient, QueryResponseRow, QueryResponseTable
from sunpy.net.scraper import Scraper
from sunpy.net.scraper_utils import get_timerange_from_exdict
from sunpy.time import TimeRange
from sunpy.util.parfive_helpers import Downloader

TIME_FORMAT = config.get("general", "time_format")

__all__ = ['QueryResponse', 'GenericClient']


class QueryResponse(QueryResponseTable):
    hide_keys = ['url']

    def time_range(self):
        """
        Returns the time-span for which records are available.
        """
        if 'Start Time' in self.colnames and 'End Time' in self.colnames:
            return TimeRange(np.min(self['Start Time']), np.max(self['End Time']))

    def response_block_properties(self):
        """
        Returns a set of class attributes on all the response blocks.
        """
        s = {a if not a.startswith('_') else None for a in dir(self[0])}
        for resp in self[1:]:
            s = s.intersection({a if not a.startswith('_') else None for a in dir(resp)})

        s.remove(None)
        return s


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

    Search uses two hooks as helper functions; these are
    :meth:`~sunpy.net.dataretriever.GenericClient.pre_search_hook` and
    :meth:`~sunpy.net.dataretriever.GenericClient.post_search_hook`.
    They help to translate the attrs for scraper before and after the search respectively.
    """
    # A string which is used to match all files and extract the desired metadata from urls correctly,
    # using ``sunpy.extern.parse.parse``.
    baseurl = None
    pattern = None
    # Set of required 'attrs' for client to handle the query.
    required = {a.Time, a.Instrument}
    # Define keywords a client needs to pass to enqueue_file
    enqueue_file_kwargs = {}

    @classmethod
    def _get_match_dict(cls, *args, **kwargs):
        """
        Constructs a dictionary using the query and registered Attrs that represents
        all possible values of the extracted metadata for files that matches the query.
        The returned dictionary is used to validate the metadata of searched files
        in :func:`~sunpy.net.scraper.Scraper._extract_files_meta`.

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
        regattrs_dict = cls.register_values()
        matchdict = {}
        for i in regattrs_dict.keys():
            attrname = i.__name__
            # only Attr values that are subclas of Simple Attr are stored as list in matchdict
            # since complex attrs like Range can't be compared with string matching.
            if issubclass(i, SimpleAttr):
                matchdict[attrname] = []
                for val, _ in regattrs_dict[i]:
                    matchdict[attrname].append(val)
        for elem in args:
            if isinstance(elem, a.Time):
                matchdict['Start Time'] = elem.start
                matchdict['End Time'] = elem.end
            elif hasattr(elem, 'value'):
                matchdict[elem.__class__.__name__] = [str(elem.value).lower()]
            elif isinstance(elem, a.Wavelength):
                matchdict['Wavelength'] = elem
            else:
                raise ValueError(
                    f"GenericClient can not add {elem.__class__.__name__} to the rowdict dictionary to pass to the Client.")
        return matchdict

    @classmethod
    def pre_search_hook(cls, *args, **kwargs):
        """
        Helper function to return the pattern and matchdict for the
        client required by :func:`~sunpy.net.dataretriever.GenericClient.search`
        before using the scraper.
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
        regattrs_dict = cls.register_values()
        optional = {k for k in regattrs_dict.keys()} - cls.required
        if not cls.check_attr_types_in_query(query, cls.required, optional):
            return False
        for key in regattrs_dict:
            all_vals = [i[0].lower() for i in regattrs_dict[key]]
            for x in query:
                if isinstance(x, key) and issubclass(key, SimpleAttr) and str(x.value).lower() not in all_vals:
                    return False
        return True

    def post_search_hook(self, exdict, matchdict):
        """
        Helper function used after :func:`~sunpy.net.dataretriever.GenericClient.search`
        which makes the extracted metadata representable in a query response table.

        Parameters
        ----------
        exdict : `dict`
            Represents metadata extracted from files.
        matchdict : `dict`
            Contains attr values accessed from ``register_values()``
            and the search query itself.

        Returns
        -------
        rowdict: `~collections.OrderedDict`
            An Ordered Dictionary which is used by `QueryResponse`
            to show results.
        """
        rowdict = OrderedDict()
        tr = get_timerange_from_exdict(exdict)
        start = tr.start
        start.format = 'iso'
        end = tr.end
        end.format = 'iso'
        rowdict['Start Time'] = start
        rowdict['End Time'] = end
        for k in matchdict:
            if k not in ('Start Time', 'End Time', 'Wavelength'):
                if k == 'Physobs':
                    # not changing case for Phsyobs
                    rowdict[k] = matchdict[k][0]
                else:
                    rowdict[k] = matchdict[k][0].upper()
        for k in exdict:
            if k not in ['year', 'month', 'day', 'hour', 'minute', 'second']:
                rowdict[k] = exdict[k]
        return rowdict

    def _get_full_filenames(self, qres, filenames, path):
        """
        Returns full pathnames for each file in the result.

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
            else:
                fname = path

            temp_dict = qres[i].response_block_map
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
        # baseurl added for backwards compatibility purposes only
        baseurl, pattern, matchdict = self.pre_search_hook(*args, **kwargs)
        scraper = Scraper(pattern=baseurl) if baseurl else Scraper(format=pattern)
        tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])
        filesmeta = scraper._extract_files_meta(tr, extractor=pattern) if baseurl else scraper._extract_files_meta(tr, matcher=matchdict)
        filesmeta = sorted(filesmeta, key=lambda k: k['url'])
        metalist = []
        for i in filesmeta:
            rowdict = self.post_search_hook(i, matchdict)
            if rowdict:
                metalist.append(rowdict)
        return QueryResponse(metalist, client=self)

    def fetch(self, qres, path=None, overwrite=False,
              progress=True, downloader=None, wait=True, **kwargs):
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
            and the existing file will be overwritten, if ``'unique'`` the filename
            will be modified to be unique.
        progress : `bool`, optional
            If `True` show a progress bar showing how many of the total files
            have been downloaded. If `False`, no progress bar will be shown.
        downloader : `parfive.Downloader`, optional
            The download manager to use.
        wait : `bool`, optional
            If `False` ``downloader.download()`` will not be called. Only has
            any effect if ``downloader`` is not `None`.

        Returns
        -------
        `parfive.Results`
        """
        if path is not None:
            path = Path(path)

        if isinstance(qres, QueryResponseRow):
            qres = qres.as_table()

        urls = []
        if len(qres):
            urls = list(qres['url'])

        filenames = [url.split('/')[-1] for url in urls]

        paths = self._get_full_filenames(qres, filenames, path)

        dl_set = True
        if not downloader:
            dl_set = False
            downloader = Downloader(progress=progress, overwrite=overwrite)

        for url, filename in zip(urls, paths):
            downloader.enqueue_file(url, filename=filename, **self.enqueue_file_kwargs)

        if dl_set and not wait:
            return

        return downloader.download()
