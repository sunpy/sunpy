# -*- coding: utf-8 -*-
"""
This module provides a wrapper around the VSO API.
"""

import os
import re
import cgi
import sys
import socket
import datetime
import warnings
import itertools
from functools import partial
from collections import defaultdict
from urllib.error import URLError, HTTPError
from urllib.request import urlopen

import zeep
from parfive import Downloader, Results
from zeep.helpers import serialize_object

import astropy.units as u
from astropy.table import QTable as Table

from sunpy import config
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient
from sunpy.net.vso import attrs
from sunpy.net.vso.attrs import _TIMEFORMAT as TIMEFORMAT
from sunpy.net.vso.attrs import _walker as walker
from sunpy.time import TimeRange, parse_time
from sunpy.util.decorators import deprecated
from sunpy.util.exceptions import SunpyUserWarning
from sunpy.util.net import get_content_disposition, slugify

from .. import _attrs as core_attrs
from .zeep_plugins import SunPyLoggingZeepPlugin

TIME_FORMAT = config.get("general", "time_format")

DEFAULT_URL_PORT = [{'url': 'http://docs.virtualsolar.org/WSDL/VSOi_rpc_literal.wsdl',
                     'port': 'nsoVSOi'},
                    {'url': 'https://sdac.virtualsolar.org/API/VSOi_rpc_literal.wsdl',
                     'port': 'sdacVSOi'}]

RANGE = re.compile(r'(\d+)(\s*-\s*(\d+))?(\s*([a-zA-Z]+))?')


# TODO: Name
class NoData(Exception):

    """ Risen for callbacks of VSOClient that are unable to supply
    information for the request. """
    pass


class _Str(str):

    """ Subclass of string that contains a meta attribute for the
    record_item associated with the file. """
    pass


# ----------------------------------------

def _parse_waverange(string):
    min_, max_, unit = RANGE.match(string).groups()[::2]
    return {
        'wave_wavemin': min_,
        'wave_wavemax': min_ if max_ is None else max_,
        'wave_waveunit': 'Angstrom' if unit is None else unit,
    }


def _parse_date(string):
    start, end = string.split(' - ')
    return {'time_start': start.strip(), 'time_end': end.strip()}


def iter_records(response):
    for prov_item in response.provideritem:
        if not hasattr(prov_item, 'record') or not prov_item.record:
            continue
        yield from prov_item.record.recorditem


def iter_errors(response):
    for prov_item in response.provideritem:
        if not hasattr(prov_item, 'record') or not prov_item.record:
            yield prov_item


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


class QueryResponse(list):
    """
    A container for VSO Records returned from VSO Searches.
    """

    def __init__(self, lst, queryresult=None, table=None):
        super().__init__(lst)
        self.queryresult = queryresult
        self.errors = []
        self.table = None

    def search(self, *query):
        """ Furtherly reduce the query response by matching it against
        another query, e.g. response.search(attrs.Instrument('aia')). """
        query = and_(*query)
        return QueryResponse(
            attrs._filter_results(query, self), self.queryresult
        )

    @classmethod
    def create(cls, queryresult):
        return cls(iter_records(queryresult), queryresult)

    def total_size(self):
        """ Total size of data in KB. May be less than the actual
        size because of inaccurate data providers. """
        # Warn about -1 values?
        return sum(record.size for record in self if record.size > 0)

    def time_range(self):
        """ Return total time-range all records span across. """
        return TimeRange(min(record.time.start for record in self if record.time.start is not None),
                         max(record.time.end for record in self if record.time.end is not None))

    def build_table(self):
        """
        Create a human readable table.

        Returns
        -------
        table : `astropy.table.QTable`
        """
        keywords = ['Start Time', 'End Time', 'Source', 'Instrument', 'Type', 'Wavelength']
        record_items = {}
        for key in keywords:
            record_items[key] = []

        def validate_time(time):
            # Handle if the time is None when coming back from VSO
            if time is None:
                return ['None']
            if record.time.start is not None:
                return [parse_time(time).strftime(TIME_FORMAT)]
            else:
                return ['N/A']

        for record in self:
            record_items['Start Time'].append(validate_time(record.time.start))
            record_items['End Time'].append(validate_time(record.time.end))
            record_items['Source'].append(str(record.source))
            record_items['Instrument'].append(str(record.instrument))
            record_items['Type'].append(str(record.extent.type)
                                        if record.extent.type is not None else ['N/A'])
            # If we have a start and end Wavelength, make a quantity
            if hasattr(record, 'wave') and record.wave.wavemin and record.wave.wavemax:
                unit = record.wave.waveunit
                # Convert this so astropy units parses it correctly
                if unit == "kev":
                    unit = "keV"
                record_items['Wavelength'].append(u.Quantity([float(record.wave.wavemin),
                                                              float(record.wave.wavemax)],
                                                             unit=unit))
            # If not save None
            else:
                record_items['Wavelength'].append(None)
        # If we have no wavelengths for the whole list, drop the col
        if all([a is None for a in record_items['Wavelength']]):
            record_items.pop('Wavelength')
            keywords.remove('Wavelength')
        else:
            # Make whole column a quantity
            try:
                with u.set_enabled_equivalencies(u.spectral()):
                    record_items['Wavelength'] = u.Quantity(record_items['Wavelength'])
            # If we have mixed units or some Nones just represent as strings
            except (u.UnitConversionError, TypeError):
                record_items['Wavelength'] = [str(a) for a in record_items['Wavelength']]

        return Table(record_items)[keywords]

    def add_error(self, exception):
        self.errors.append(exception)

    def response_block_properties(self):
        """
        Returns a set of class attributes on all the response blocks.

        Returns
        -------
        s : `set`
            List of strings, containing attribute names in the response blocks.
        """
        s = {a if not a.startswith('_') else None for a in dir(self[0])}
        for resp in self[1:]:
            s = s.intersection({a if not a.startswith('_') else None for a in dir(resp)})

        s.remove(None)
        return s

    def __str__(self):
        """Print out human-readable summary of records retrieved"""
        return str(self.build_table())

    def __repr__(self):
        """Print out human-readable summary of records retrieved"""
        return repr(self.build_table())

    def _repr_html_(self):
        return self.build_table()._repr_html_()


class DownloadFailed(Exception):
    pass


class MissingInformation(Exception):
    pass


class UnknownMethod(Exception):
    pass


class MultipleChoices(Exception):
    pass


class UnknownVersion(Exception):
    pass


class UnknownStatus(Exception):
    pass


class VSOClient(BaseClient):
    """
    VSO Client

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

    def make(self, atype, **kwargs):
        """
        Create a new SOAP object.
        """
        obj = self.api.get_type(f"VSO:{atype}")
        return obj(**kwargs)

    def search(self, *query):
        """ Query data from the VSO with the new API. Takes a variable number
        of attributes as parameter, which are chained together using AND.

        The new query language allows complex queries to be easily formed.

        Examples
        --------
        Query all data from eit or aia between 2010-01-01T00:00 and
        2010-01-01T01:00.

        >>> from datetime import datetime
        >>> from sunpy.net import vso, attrs as a
        >>> client = vso.VSOClient()  # doctest: +REMOTE_DATA
        >>> client.search(
        ...    a.Time(datetime(2010, 1, 1), datetime(2010, 1, 1, 1)),
        ...    a.Instrument('eit') | a.Instrument('aia'))   # doctest:  +REMOTE_DATA
        <QTable length=5>
           Start Time [1]       End Time [1]    Source ...   Type   Wavelength [2]
                                                       ...             Angstrom
               str19               str19         str4  ...   str8      float64
        ------------------- ------------------- ------ ... -------- --------------
        2010-01-01 00:00:08 2010-01-01 00:00:20   SOHO ... FULLDISK 195.0 .. 195.0
        2010-01-01 00:12:08 2010-01-01 00:12:20   SOHO ... FULLDISK 195.0 .. 195.0
        2010-01-01 00:24:10 2010-01-01 00:24:22   SOHO ... FULLDISK 195.0 .. 195.0
        2010-01-01 00:36:08 2010-01-01 00:36:20   SOHO ... FULLDISK 195.0 .. 195.0
        2010-01-01 00:48:09 2010-01-01 00:48:21   SOHO ... FULLDISK 195.0 .. 195.0

        Returns
        -------
        out : :py:class:`QueryResult` (enhanced list)
            Matched items. Return value is of same type as the one of
            :py:meth:`VSOClient.search`.
        """
        query = and_(*query)
        QueryRequest = self.api.get_type('VSO:QueryRequest')
        VSOQueryResponse = self.api.get_type('VSO:QueryResponse')
        responses = []
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
                response = QueryResponse.create(self.merge(responses))
                response.add_error(ex)

        return QueryResponse.create(self.merge(responses))

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
    def mk_filename(pattern, queryresponse, resp, url):
        """
        Generate the best possible (or least-worse) filename for a VSO download.

        * Use the ``content-disposition`` header.
        * Use `fileid` to generate a file name if content-disposition fails
        * If everything else fails use the last segment of the URL and hope.
        """
        name = None
        if resp:
            cdheader = resp.headers.get("Content-Disposition", None)
            if cdheader:
                value, params = cgi.parse_header(cdheader)
                name = params.get('filename', "")
                # Work around https://github.com/sunpy/sunpy/issues/3372
                if name.count('"') >= 2:
                    name = name.split('"')[1]

        if name is None:
            # Advice from the VSO is to fallback to providerid + fileid
            # As it's possible multiple providers give the same fileid.
            # However, I haven't implemented this yet as it would be a breaking
            # change to the filenames we expect.

            # I don't know if we still need this bytes check in Python 3 only
            # land, but I don't dare remove it.
            if isinstance(queryresponse.fileid, bytes):
                fileid = queryresponse.fileid.decode("ascii", "ignore")
            else:
                fileid = queryresponse.fileid

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

        fname = pattern.format(file=name, **serialize_object(queryresponse))

        return fname

    @deprecated("1.0", alternative="sunpy.net.Fido")
    def query_legacy(self, tstart=None, tend=None, **kwargs):
        """
        Query data from the VSO mocking the IDL API as close as possible.
        Either tstart and tend or date_start and date_end or date have
        to be supplied.

        Parameters
        ----------
        tstart : datetime.datetime
            Start of the time-range in which records are searched.
        tend : datetime.datetime
            Start of the time-range in which records are searched.
        date : str
            (start date) - (end date)
        start_date : datetime
            the start date
        end_date : datetime
            the end date
        wave : str
            (min) - (max) (unit)
        min_wave : str
            minimum spectral range
        max_wave : str
            maximum spectral range
        unit_wave : str
            spectral range units (Angstrom, GHz, keV)
        extent : str
            VSO 'extent type' ... (FULLDISK, CORONA, LIMB, etc)
        physobj : str
            VSO 'physical observable'
        provider : str
            VSO ID for the data provider (SDAC, NSO, SHA, MSU, etc)
        source : str
            spacecraft or observatory (SOHO, YOHKOH, BBSO, etc)
            synonyms : spacecraft, observatory
        instrument : str
            instrument ID (EIT, SXI-0, SXT, etc)
            synonyms : telescope, inst
        detector : str
            detector ID (C3, EUVI, COR2, etc.)
        layout : str
            layout of the data (image, spectrum, time_series, etc.)
        level : str
            level of the data product (numeric range, see below)
        pixels : str
            number of pixels (numeric range, see below)
        resolution : str
            effective resolution (1 = full, 0.5 = 2x2 binned, etc)
            numeric range, see below.
        pscale : str
            pixel scale, in arcseconds (numeric range, see below)
        near_time : datetime
            return record closest to the time.  See below.
        sample : int
            attempt to return only one record per SAMPLE seconds.  See below.

        Numeric Ranges:

            - May be entered as a string or any numeric type for equality matching
            - May be a string of the format '(min) - (max)' for range matching
            - May be a string of the form '(operator) (number)' where operator
              is one of: lt gt le ge < > <= >=


        Examples
        --------
        Query all data from eit between 2010-01-01T00:00 and
        2010-01-01T01:00.

        >>> from datetime import datetime
        >>> from sunpy.net import vso
        >>> client = vso.VSOClient()  # doctest: +SKIP
        >>> qr = client.query_legacy(datetime(2010, 1, 1),
        ...                          datetime(2010, 1, 1, 1),
        ...                          instrument='eit')  # doctest: +SKIP

        Returns
        -------
        out : :py:class:`QueryResult` (enhanced list)
            Matched items. Return value is of same type as the one of
            :py:class:`VSOClient.search`.
        """
        def sdk(key): return partial(lambda key, value: {key: value}, key)
        ALIASES = {
            'wave_min': sdk('wave_wavemin'),
            'wave_max': sdk('wave_wavemax'),
            'wave_type': sdk('wave_wavetype'),
            'wave_unit': sdk('wave_waveunit'),
            'min_wave': sdk('wave_wavemin'),
            'max_wave': sdk('wave_wavemax'),
            'type_wave': sdk('wave_wavetype'),
            'unit_wave': sdk('wave_waveunit'),
            'wave': _parse_waverange,
            'inst': sdk('instrument'),
            'telescope': sdk('instrument'),
            'spacecraft': sdk('source'),
            'observatory': sdk('source'),
            'start_date': sdk('time_start'),
            'end_date': sdk('time_end'),
            'start': sdk('time_start'),
            'end': sdk('time_end'),
            'near_time': sdk('time_near'),
            'date': _parse_date,
            'layout': sdk('datatype'),
        }
        if tstart is not None:
            kwargs.update({'time_start': tstart})
        if tend is not None:
            kwargs.update({'time_end': tend})

        QueryRequest = self.api.get_type('VSO:QueryRequest')
        VSOQueryResponse = self.api.get_type('VSO:QueryResponse')
        block = self.api.get_type('VSO:QueryRequestBlock')()

        for key, value in kwargs.items():
            for k, v in ALIASES.get(key, sdk(key))(value).items():
                if k.startswith('time'):
                    v = parse_time(v).strftime(TIMEFORMAT)
                attr = k.split('_')
                lst = attr[-1]
                rest = attr[:-1]

                for elem in rest:
                    try:
                        if block[elem] is None:
                            block[elem] = {}
                        block = block[elem]
                    except KeyError:
                        raise ValueError(
                            f"Unexpected argument {key!s}.")
                if lst in block and block[lst]:
                    raise ValueError(
                        f"Got multiple values for {k!s}.")
                block[lst] = v

        return QueryResponse.create(VSOQueryResponse(
            self.api.service.Query(QueryRequest(block=block))))

    @deprecated("1.0")
    def latest(self):
        """ Return newest record (limited to last week). """
        from datetime import datetime, timedelta
        return self.query_legacy(
            datetime.utcnow() - timedelta(7),
            datetime.utcnow(),
            time_near=datetime.utcnow()
        )

    def fetch(self, query_response, path=None, methods=None, site=None,
              progress=True, overwrite=False, downloader=None, wait=True):
        """
        Download data specified in the query_response.

        Parameters
        ----------
        query_response : sunpy.net.vso.QueryResponse
            QueryResponse containing the items to be downloaded.

        path : str
            Specify where the data is to be downloaded. Can refer to arbitrary
            fields of the QueryResponseItem (instrument, source, time, ...) via
            string formatting, moreover the file-name of the file downloaded can
            be referred to as file, e.g.
            "{source}/{instrument}/{time.start}/{file}".

        methods : {list of str}
            Download methods, defaults to URL-FILE_Rice then URL-FILE.
            Methods are a concatenation of one PREFIX followed by any number of
            SUFFIXES i.e. `PREFIX-SUFFIX_SUFFIX2_SUFFIX3`.
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
            and the existing file will be overwritten, if `'unique'` the filename
            will be modified to be unique.

        downloader : `parfive.Downloader`, optional
            The download manager to use.

        wait : `bool`, optional
           If `False` ``downloader.download()`` will not be called. Only has
           any effect if `downloader` is not `None`.

        Returns
        -------
        out : `parfive.Results`
            Object that supplies a list of filenames and any errors.

        Examples
        --------
        >>> files = fetch(qr) # doctest:+SKIP
        """
        if path is None:
            path = os.path.join(config.get('downloads', 'download_dir'),
                                '{file}')
        elif isinstance(path, str) and '{file}' not in path:
            path = os.path.join(path, '{file}')
        path = os.path.expanduser(path)

        dl_set = True
        if not downloader:
            dl_set = False
            downloader = Downloader(progress=progress)

        fileids = VSOClient.by_fileid(query_response)
        if not fileids:
            return downloader.download() if wait else Results()
        # Adding the site parameter to the info
        info = {}
        if site is not None:
            info['site'] = site

        VSOGetDataResponse = self.api.get_type("VSO:VSOGetDataResponse")

        data_request = self.make_getdatarequest(query_response, methods, info)
        data_response = VSOGetDataResponse(self.api.service.GetData(data_request))

        err_results = self.download_all(data_response, methods, downloader, path, fileids)

        if dl_set and not wait:
            return err_results

        results = downloader.download()
        results += err_results
        results._errors += err_results.errors
        return results

    @staticmethod
    def link(query_response, maps):
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
            # pylint: disable=W0201
            item.meta = record_item
            ret.append(item)
        return ret

    def make_getdatarequest(self, response, methods=None, info=None):
        """ Make datarequest with methods from response. """
        if methods is None:
            methods = self.method_order + ['URL']

        return self.create_getdatarequest(
            {k: [x.fileid for x in v]
                 for k, v in self.by_provider(response).items()},
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

    # pylint: disable=R0913,R0912
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
            # pylint: disable=W0631
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
    def by_provider(response):
        """
        Returns a dictionary of provider
        corresponding to records in the response.
        """

        map_ = defaultdict(list)
        for record in response:
            map_[record.provider].append(record)
        return map_

    @staticmethod
    def by_fileid(response):
        """
        Returns a dictionary of fileids
        corresponding to records in the response.
        """

        return {
            record.fileid: record for record in response
        }

    # pylint: disable=W0613
    def multiple_choices(self, choices, response):
        """ Override to pick between multiple download choices. """
        for elem in self.method_order:
            if elem in choices:
                return [elem]
        raise NoData

    # pylint: disable=W0613
    def missing_information(self, info, field):
        """ Override to provide missing information. """
        raise NoData

    # pylint: disable=W0613
    def unknown_method(self, response):
        """ Override to pick a new method if the current one is unknown. """
        raise NoData

    @classmethod
    def _can_handle_query(cls, *query):
        return all([x.__class__.__name__ in core_attrs.__all__ + attrs.__all__ for x in query])
