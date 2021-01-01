"""
This module provides a wrapper around the VSO API.
"""

import os
import cgi
import copy
import json
import socket
import inspect
import datetime
import warnings
import itertools
from pathlib import Path
from functools import partial
from urllib.error import URLError, HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import zeep

from sunpy import config, log
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseRow
from sunpy.net.vso import attrs
from sunpy.net.vso.attrs import _walker as walker
from sunpy.util.decorators import deprecated
from sunpy.util.exceptions import SunpyDeprecationWarning, SunpyUserWarning
from sunpy.util.net import slugify
from sunpy.util.parfive_helpers import Downloader, Results
from .. import _attrs as core_attrs
from .exceptions import (
    DownloadFailed,
    MissingInformation,
    MultipleChoices,
    NoData,
    UnknownMethod,
    UnknownStatus,
    UnknownVersion,
)
from .legacy_response import QueryResponse
from .table_response import VSOQueryResponseTable
from .zeep_plugins import SunPyLoggingZeepPlugin

DEFAULT_URL_PORT = [{'url': 'http://docs.virtualsolar.org/WSDL/VSOi_rpc_literal.wsdl',
                     'port': 'nsoVSOi'},
                    {'url': 'https://sdac.virtualsolar.org/API/VSOi_rpc_literal.wsdl',
                     'port': 'sdacVSOi'}]


class _Str(str):
    """ Subclass of string that contains a meta attribute for the
    record_item associated with the file. """
    meta = None


# ----------------------------------------
def check_connection(url):
    try:
        return urlopen(url).getcode() == 200
    except (socket.error, socket.timeout, HTTPError, URLError) as e:
        warnings.warn(f"Connection to {url} failed with error {e}. Retrying with different url and port.",
                      SunpyUserWarning)
        return None


def get_online_vso_url():
    """
    Return the first VSO url and port combination that is online.
    """
    for mirror in DEFAULT_URL_PORT:
        if check_connection(mirror['url']):
            return mirror


def build_client(url=None, port_name=None, **kwargs):
    """
    Construct a `zeep.Client` object to connect to VSO.

    Parameters
    ----------
    url : `str`
        The URL to connect to.

    port_name : `str`
        The "port" to use.

    kwargs : `dict`
        All extra keyword arguments are passed to `zeep.Client`.

    Returns
    -------

    `zeep.Client`
    """
    if url is None and port_name is None:
        mirror = get_online_vso_url()
        if mirror is None:
            raise ConnectionError("No online VSO mirrors could be found.")
        url = mirror['url']
        port_name = mirror['port']
    elif url and port_name:
        if not check_connection(url):
            raise ConnectionError(f"Can't connect to url {url}")
    else:
        raise ValueError("Both url and port_name must be specified if either is.")

    if "plugins" not in kwargs:
        kwargs["plugins"] = [SunPyLoggingZeepPlugin()]

    client = zeep.Client(url, port_name=port_name, **kwargs)
    client.set_ns_prefix('VSO', 'http://virtualsolar.org/VSO/VSOi')
    return client


class VSOClient(BaseClient):
    """
    Provides access to query and download from Virtual Solar Observatory (VSO).

    Parameters
    ----------
    url : `str`, optional
        The VSO url to use. If not specified will use the first online known URL.

    port : `str`, optional
        The VSO port name to use. If not specified will use the first online known URL.

    api : `zeep.Client`, optional
        The `zeep.Client` instance to use for interacting with the VSO. If not
        specified one will be created.
    """
    method_order = [
        'URL-FILE_Rice', 'URL-FILE', 'URL-packaged', 'URL-TAR_GZ', 'URL-ZIP', 'URL-TAR',
    ]

    def __init__(self, url=None, port=None, api=None):
        if not isinstance(api, zeep.Client):
            api = build_client(url, port)
            if api is None:
                raise ConnectionError("Cannot find an online VSO mirror.")
        self.api = api

    def __deepcopy__(self, memo):
        """
        Copy the client but don't copy the API object.
        """
        memo[id(self.api)] = self.api
        deepcopy_method = self.__deepcopy__
        self.__deepcopy__ = None
        cp = copy.deepcopy(self, memo)
        self.__deepcopy__ = deepcopy_method
        cp.__deepcopy__ = deepcopy_method
        return cp

    def make(self, atype, **kwargs):
        """
        Create a new SOAP object.
        """
        obj = self.api.get_type(f"VSO:{atype}")
        return obj(**kwargs)

    def search(self, *query, response_format=None):
        """
        Query data from the VSO with the new API. Takes a variable number
        of attributes as parameter, which are chained together using AND.

        Parameters
        ----------
        response_format: {"legacy", "table"}
            The response format from the search, this can be either
            ``"legacy"`` to return a list-like object of the zeep responses, or
            ``"table"`` to return the responses in a subclass of
            `~astropy.table.QTable`.

        Examples
        --------
        Query all data from eit or aia between 2010-01-01T00:00 and
        2010-01-01T01:00.

        >>> from datetime import datetime
        >>> from sunpy.net import vso, attrs as a
        >>> client = vso.VSOClient()  # doctest: +REMOTE_DATA
        >>> client.search(
        ...    a.Time(datetime(2010, 1, 1), datetime(2010, 1, 1, 1)),
        ...    a.Instrument.eit | a.Instrument.aia,
        ...    response_format="table")   # doctest:  +REMOTE_DATA
        <sunpy.net.vso.table_response.VSOQueryResponseTable object at ...>
            Start Time               End Time        Source ... Extent Type   Size
                                                            ...              Mibyte
        ----------------------- ----------------------- ------ ... ----------- -------
        2010-01-01 00:00:08.000 2010-01-01 00:00:20.000   SOHO ...    FULLDISK 2.01074
        2010-01-01 00:12:08.000 2010-01-01 00:12:20.000   SOHO ...    FULLDISK 2.01074
        2010-01-01 00:24:10.000 2010-01-01 00:24:22.000   SOHO ...    FULLDISK 2.01074
        2010-01-01 00:36:08.000 2010-01-01 00:36:20.000   SOHO ...    FULLDISK 2.01074
        2010-01-01 00:48:09.000 2010-01-01 00:48:21.000   SOHO ...    FULLDISK 2.01074


        Returns
        -------
        out : `QueryResult`
            Matched items. Return value is of same type as the one of
            :meth:`VSOClient.search`.
        """
        if response_format is None:
            response_format = "legacy"
            warnings.warn("The default response format from the VSO client will "
                          "be changing to 'table' in version 3.1. "
                          "To remove this warning set response_format='legacy' "
                          "to maintain the old behaviour or response_format='table'"
                          " to use the new behaviour.",
                          SunpyDeprecationWarning,
                          stacklevel=2)
        query = and_(*query)
        QueryRequest = self.api.get_type('VSO:QueryRequest')
        VSOQueryResponse = self.api.get_type('VSO:QueryResponse')
        responses = []
        exceptions = []
        for block in walker.create(query, self.api):
            try:
                query_response = self.api.service.Query(
                    QueryRequest(block=block)
                )
                for resp in query_response:
                    if resp["error"]:
                        warnings.warn(resp["error"], SunpyUserWarning)
                responses.append(
                    VSOQueryResponse(query_response)
                )
            except Exception as ex:
                exceptions.append(ex)

        responses = self.merge(responses)
        if response_format == "legacy":
            response = QueryResponse.create(responses)
        else:
            response = VSOQueryResponseTable.from_zeep_response(responses, client=self)

        for ex in exceptions:
            response.add_error(ex)

        return response

    def merge(self, queryresponses):
        """ Merge responses into one. """
        if len(queryresponses) == 1:
            return queryresponses[0]

        fileids = set()
        providers = {}

        for queryresponse in queryresponses:
            for provideritem in queryresponse.provideritem:
                provider = provideritem.provider
                if not hasattr(provideritem, 'record'):
                    continue
                if not hasattr(provideritem.record, 'recorditem'):
                    continue
                if provideritem.provider not in providers:
                    providers[provider] = provideritem
                    fileids |= {
                        record_item.fileid
                        for record_item in provideritem.record.recorditem
                    }
                else:
                    for record_item in provideritem.record.recorditem:
                        if record_item.fileid not in fileids:
                            fileids.add(record_item.fileid)
                            providers[provider].record.recorditem.append(
                                record_item
                            )
                            providers[provider].no_of_records_found += 1
                            providers[provider].no_of_records_returned += 1
        return self.make('QueryResponse',
                         provideritem=list(providers.values()))

    @staticmethod
    def mk_filename(pattern, queryresponserow, resp, url):
        """
        Generate the best possible (or least-worse) filename for a VSO download.

        * Use the ``content-disposition`` header.
        * Use ``fileid`` to generate a file name if content-disposition fails
        * If everything else fails use the last segment of the URL and hope.
        """
        name = None
        if resp:
            cdheader = resp.headers.get("Content-Disposition", None)
            if cdheader:
                _, params = cgi.parse_header(cdheader)
                name = params.get('filename', "")
                # Work around https://github.com/sunpy/sunpy/issues/3372
                if name.count('"') >= 2:
                    name = name.split('"')[1]

        if name is None:
            # Advice from the VSO is to fallback to providerid + fileid for a filename
            # As it's possible multiple providers give the same fileid.
            # However, I haven't implemented this yet as it would be a breaking
            # change to the filenames we expect.
            fileid = queryresponserow['fileid']

            # Some providers make fileid a path
            # Some also don't specify a file extension, but not a lot we can do
            # about that.
            name = fileid.split("/")[-1]

        # If somehow we have got this far with an empty string, fallback to url segment
        if not name:
            name = url.split('/')[-1]

        # Remove any not-filename appropriate characters
        name = slugify(name)

        # If absolutely everything else fails make a filename based on download time
        if not name:
            name = f"vso_file_{datetime.datetime.now().strftime('%Y%m%d%H%M%S%f')}"

        fname = pattern.format(file=name,
                               **queryresponserow.response_block_map)

        return fname

    def fetch(self, query_response, path=None, methods=None, site=None,
              progress=True, overwrite=False, downloader=None, wait=True):
        """
        Download data specified in the query_response.

        Parameters
        ----------
        query_response : sunpy.net.vso.VSOQueryResponseTable
            QueryResponse containing the items to be downloaded.

        path : str
            Specify where the
        data is to be downloaded. Can refer to arbitrary
            fields of the QueryResponseItem (instrument, source, time, ...) via
            string formatting, moreover the file-name of the file downloaded can
            be referred to as file, e.g.
            "{source}/{instrument}/{time.start}/{file}".

        methods : {list of str}
            Download methods, defaults to URL-FILE_Rice then URL-FILE.
            Methods are a concatenation of one PREFIX followed by any number of
            SUFFIXES i.e. ``PREFIX-SUFFIX_SUFFIX2_SUFFIX3``.
            The full list of
            `PREFIXES <https://sdac.virtualsolar.org/cgi/show_details?keyword=METHOD_PREFIX>`_
            and `SUFFIXES <https://sdac.virtualsolar.org/cgi/show_details?keyword=METHOD_SUFFIX>`_
            are listed on the VSO site.

        site : str
            There are a number of caching mirrors for SDO and other
            instruments, some available ones are listed below.

            =============== ========================================================
            NSO             National Solar Observatory, Tucson (US)
            SAO  (aka CFA)  Smithonian Astronomical Observatory, Harvard U. (US)
            SDAC (aka GSFC) Solar Data Analysis Center, NASA/GSFC (US)
            ROB             Royal Observatory of Belgium (Belgium)
            MPS             Max Planck Institute for Solar System Research (Germany)
            UCLan           University of Central Lancashire (UK)
            IAS             Institut Aeronautique et Spatial (France)
            KIS             Kiepenheuer-Institut fur Sonnenphysik Germany)
            NMSU            New Mexico State University (US)
            =============== ========================================================

        progress : `bool`, optional
            If `True` show a progress bar showing how many of the total files
            have been downloaded. If `False`, no progress bars will be shown at all.

        overwrite : `bool` or `str`, optional
            Determine how to handle downloading if a file already exists with the
            same name. If `False` the file download will be skipped and the path
            returned to the existing file, if `True` the file will be downloaded
            and the existing file will be overwritten, if ``'unique'`` the filename
            will be modified to be unique.

        downloader : `parfive.Downloader`, optional
            The download manager to use.

        wait : `bool`, optional
           If `False` ``downloader.download()`` will not be called. Only has
           any effect if ``downloader`` is not `None`.

        Returns
        -------
        out : `parfive.Results`
            Object that supplies a list of filenames and any errors.

        Examples
        --------
        >>> files = fetch(qr) # doctest:+SKIP
        """
        if path is None:
            path = Path(config.get('downloads', 'download_dir')) / '{file}'
        elif isinstance(path, (str, os.PathLike)) and '{file}' not in str(path):
            path = Path(path) / '{file}'
        else:
            path = Path(path)
        path = path.expanduser()

        dl_set = True
        if not downloader:
            dl_set = False
            downloader = Downloader(progress=progress, overwrite=overwrite)

        if isinstance(query_response, (QueryResponse, list)):
            query_response = VSOQueryResponseTable.from_zeep_response(query_response,
                                                                      client=self,
                                                                      _sort=False)
        if isinstance(query_response, QueryResponseRow):
            query_response = query_response.as_table()

        if not len(query_response):
            return downloader.download() if wait else Results()

        # Adding the site parameter to the info
        info = {}
        if site is not None:
            info['site'] = site

        VSOGetDataResponse = self.api.get_type("VSO:VSOGetDataResponse")

        data_request = self.make_getdatarequest(query_response, methods, info)
        data_response = VSOGetDataResponse(self.api.service.GetData(data_request))

        err_results = self.download_all(data_response,
                                        methods,
                                        downloader,
                                        str(path),
                                        self.by_fileid(query_response))

        if dl_set and not wait:
            return err_results

        results = downloader.download()
        results += err_results
        results._errors += err_results.errors
        return results

    @deprecated("2.1", "This functionality is deprecated as it is replaced by better search support.")
    @staticmethod
    def link(query_response, maps):  # pragma: no cover
        """ Return list of paths with records associated with them in
        the meta attribute. """
        if not maps:
            return []
        ret = []

        for record_item in query_response:
            try:
                item = _Str(maps[record_item.fileid]['path'])
            except KeyError:
                continue
            item.meta = record_item
            ret.append(item)
        return ret

    def make_getdatarequest(self, response, methods=None, info=None):
        """ Make datarequest with methods from response. """
        if methods is None:
            methods = self.method_order + ['URL']

        return self.create_getdatarequest(
            {g[0]['Provider']: list(g['fileid']) for g in response.group_by('Provider').groups},
            methods, info
        )

    def create_getdatarequest(self, maps, methods, info=None):
        """ Create datarequest from maps mapping data provider to
        fileids and methods, """
        if info is None:
            info = {}

        if 'email' not in info:
            info['email'] = 'sunpy'

        # For the JSOC provider we need to make a DataRequestItem for each
        # series, not just one for the whole provider.

        # Remove JSOC provider items from the map
        jsoc = maps.pop('JSOC', [])
        # Make DRIs for everything that's not JSOC one per provider
        dris = [self.make('DataRequestItem', provider=k, fileiditem={'fileid': v})
                for k, v in maps.items()]

        def series_func(x):
            """ Extract the series from the fileid. """
            return x.split(':')[0]

        # Sort the JSOC fileids by series
        # This is a precursor to groupby as recommended by the groupby docs
        series_sorted = sorted(jsoc, key=series_func)
        # Iterate over the series and make a DRI for each.
        # groupby creates an iterator based on a key function, in this case
        # based on the series (the part before the first ':')
        for series, fileids in itertools.groupby(series_sorted, key=series_func):
            dris.append(self.make('DataRequestItem',
                                  provider='JSOC',
                                  fileiditem={'fileid': list(fileids)}))

        request = {'method': {'methodtype': methods},
                   'info': info,
                   'datacontainer': {'datarequestitem': dris}
                   }

        return self.make('VSOGetDataRequest', request=request)

    def download_all(self, response, methods, downloader, path, qr, info=None):
        results = Results()
        GET_VERSION = [
            ('0.8', (5, 8)),
            ('0.7', (1, 4)),
            ('0.6', (0, 3)),
        ]

        for dresponse in response.getdataresponseitem:
            for version, (from_, to) in GET_VERSION:
                if getattr(dresponse, version, '0.6') >= version:
                    break
            else:
                results.add_error('', UnknownVersion(dresponse))
                continue

            # If from_ and to are uninitialized, the else block of the loop
            # continues the outer loop and thus this code is never reached.
            code = (
                dresponse.status[from_:to]
                if getattr(dresponse, 'status', None) else '200'
            )
            if code == '200':
                for dataitem in dresponse.getdataitem.dataitem:

                    try:
                        self.download(
                            dresponse.method.methodtype[0],
                            dataitem.url,
                            downloader,
                            path,
                            qr[dataitem.fileiditem.fileid[0]]
                        )
                    except NoData:
                        results.add_error('', '', DownloadFailed(dresponse))
                        continue

            elif code == '300' or code == '412' or code == '405':
                if code == '300':
                    try:
                        methods = self.multiple_choices(
                            dresponse.method.methodtype, dresponse
                        )
                    except NoData:
                        results.add_error('', '', MultipleChoices(dresponse))
                        continue
                elif code == '412':
                    try:
                        info = self.missing_information(
                            info, dresponse.info
                        )
                    except NoData:
                        results.add_error('', '', MissingInformation(dresponse))
                        continue
                elif code == '405':
                    try:
                        methods = self.unknown_method(dresponse)
                    except NoData:
                        results.add_error('', '', UnknownMethod(dresponse))
                        continue

                files = []
                for dataitem in dresponse.getdataitem.dataitem:
                    files.extend(dataitem.fileiditem.fileid)

                request = self.create_getdatarequest(
                    {dresponse.provider: files}, methods, info
                )

                self.download_all(
                    self.api.service.GetData(request), methods, downloader, path,
                    qr, info
                )
            else:
                results.add_error('', '', UnknownStatus(dresponse))

        return results

    def download(self, method, url, downloader, *args):
        """ Enqueue a file to be downloaded, extra args are passed to ``mk_filename``"""
        if method.startswith('URL'):
            return downloader.enqueue_file(url, filename=partial(self.mk_filename, *args))

        raise NoData

    @staticmethod
    def by_fileid(response):
        """
        Returns a dictionary of fileids
        corresponding to records in the response.
        """
        return {
            record['fileid']: record for record in response
        }

    def multiple_choices(self, choices, response):
        """ Override to pick between multiple download choices. """
        for elem in self.method_order:
            if elem in choices:
                return [elem]
        raise NoData

    def missing_information(self, info, field):
        """ Override to provide missing information. """
        raise NoData

    def unknown_method(self, response):
        """ Override to pick a new method if the current one is unknown. """
        raise NoData

    @classmethod
    def _can_handle_query(cls, *query):
        required = {core_attrs.Time}
        # Get all classes in core_attrs and attrs
        optional = {value for (name, value) in inspect.getmembers(core_attrs) if
                    name in core_attrs.__all__}
        optional.update(value for (name, value) in inspect.getmembers(attrs) if
                        name in attrs.__all__)
        return cls.check_attr_types_in_query(query, required, optional)

    @classmethod
    def _attrs_module(cls):
        return 'vso', 'sunpy.net.vso.attrs'

    def __del__(self):
        """
        Attempt to close the connection, but if it fails, continue.
        """
        try:
            self.api.transport.session.close()
        except Exception as e:
            log.debug(f"Failed to close VSO API connection with: {e}")

    @classmethod
    def register_values(cls):
        # We always use the local file for now.
        return cls.load_vso_values()

    @staticmethod
    def load_vso_values():
        """
        We take this list and register all the keywords as corresponding Attrs.

        Returns
        -------
        dict
            The constructed Attrs dictionary ready to be passed into Attr registry.
        """
        from sunpy.net import attrs as a

        here = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(here, 'data', 'attrs.json'), 'r') as attrs_file:
            keyword_info = json.load(attrs_file)

        # Now to traverse the saved dict and give them attr keys.
        attrs = {}
        for key, value in keyword_info.items():
            attr = getattr(a, key.capitalize(), None)
            if attr is None:
                attr = getattr(a.vso, key.capitalize())
            attrs[attr] = value
        return attrs

    @staticmethod
    def create_parse_vso_values():
        """
        Makes a network call to the VSO API that returns what keywords they support.
        We take this list and register all the keywords as corresponding Attrs.
        """
        here = os.path.dirname(os.path.realpath(__file__))

        # Keywords we are after
        keywords = ["+detector", "+instrument", "+source", "+provider", "+physobs", "+level"]
        # Construct and format the request
        keyword_info = {}
        url = "https://vso1.nascom.nasa.gov/cgi-bin/registry_json.cgi"
        headers = {"Content-Type": "application/x-www-form-urlencoded"}
        for keyword in keywords:
            data = urlencode({'fields': f"['{keyword}']".replace("'", '"')}).encode('ascii')
            req = Request(url=url, data=data, headers=headers)
            response = urlopen(req)
            keyword_info[keyword.replace("+", "")] = json.loads(response.read())

        # Now to traverse the return and create attrs out of them.
        attrs = {}
        for key, value in keyword_info.items():
            attrs[key] = []
            for item in value:
                if item:
                    if key == "level":
                        attrs[key].append((str(item[key]), str(item[key])))
                    else:
                        attrs[key].append((str(item[key]), str(item[key+"_long"])))

        with open(os.path.join(here, 'data', 'attrs.json'), 'w') as attrs_file:
            json.dump(attrs, attrs_file, indent=2)
