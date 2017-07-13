# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import

import os
import time
import warnings

import requests
import numpy as np
import astropy.units as u
import astropy.time
import astropy.table
from astropy.utils.misc import isiterable
import drms

from sunpy import config
from sunpy.time import parse_time, TimeRange
from sunpy.net.download import Downloader, Results
from sunpy.net.attr import and_
from sunpy.net.jsoc import attrs
from sunpy.net.jsoc.attrs import walker, Keys, PrimeKeys
from sunpy.net.vso.attrs import _VSOSimpleAttr
from sunpy.extern.six.moves import urllib
from sunpy.extern import six
from sunpy.util.metadata import MetaDict
from sunpy.util import deprecated

__all__ = ['JSOCClient', 'JSOCResponse']


class JSOCResponse(object):
    def __init__(self, table=None):
        """
        table : `astropy.table.Table`
        """

        self.table = table
        self.query_args = None
        self.requestIDs = None

    def __str__(self):
        return str(self.table)

    def __repr__(self):
        return repr(self.table)

    def _repr_html_(self):
        return self.table._repr_html_()

    def __len__(self):
        if self.table is None:
            return 0
        else:
            return len(self.table)

    def append(self, table):
        if self.table is None:
            self.table = table
        else:
            self.table = astropy.table.vstack([self.table, table])


class JSOCClient(object):
    """
    This is a Client to the JSOC Data Export service.

    It exposes a similar API to the VSO client, although the underlying model
    is more complex. The JSOC stages data before you can download it, so a JSOC
    query is a three stage process, first you query the JSOC for records,
    a table of these records is returned. Then you can request these records to
    be staged for download and then you can download them.
    The last two stages of this process are bundled together into the `get()`
    method, but they can be separated if you are performing a large or complex
    query.

    .. warning::
        JSOC now requires you to register your email address before requesting
        data. See this site: http://jsoc.stanford.edu/ajax/register_email.html

    Notes
    -----
    This Client mocks input to this site: http://jsoc.stanford.edu/ajax/exportdata.html
    Therefore that is a good resource if things are mis-behaving.
    The full list of 'series' is available through this site: http://jsoc.stanford.edu/

    You can build more complex queries by specifying parameters to POST to JSOC via keyword
    arguments. You can generate these kwargs using the Export Data page at JSOC.

    JSOC now requires a validated email address, you can pass in your validated email address
    using the `~sunpy.net.jsoc.attrs.Notify` attribute. You have to register your email address
    with JSOC http://jsoc.stanford.edu/ajax/register_email.html.


    Examples
    --------

    *Example 1*

    Query JSOC for some HMI data at 45 second cadence:

    >>> from sunpy.net import jsoc
    >>> from sunpy.net import attrs as a
    >>> client = jsoc.JSOCClient()
    >>> response = client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                         a.jsoc.Series('hmi.m_45s'), a.jsoc.Notify("sunpy@sunpy.org"))

    the response object holds the records that your query will return:

    >>> print(response)   # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
            DATE         TELESCOP  INSTRUME  ... WAVELNTH     WAVEUNIT
    -------------------- -------- ---------- ... -------- ---------------
    2014-01-05T17:44:53Z  SDO/HMI HMI_FRONT2 ...   6173.0 Invalid KeyLink
    2014-01-05T17:46:02Z  SDO/HMI HMI_FRONT2 ...   6173.0 Invalid KeyLink
    2014-01-05T17:47:11Z  SDO/HMI HMI_FRONT2 ...   6173.0 Invalid KeyLink
    2014-01-05T17:48:18Z  SDO/HMI HMI_FRONT2 ...   6173.0 Invalid KeyLink
                     ...      ...        ... ...      ...             ...
    2014-01-05T17:42:33Z  SDO/HMI HMI_FRONT2 ...   6173.0 Invalid KeyLink
    2014-01-05T17:43:41Z  SDO/HMI HMI_FRONT2 ...   6173.0 Invalid KeyLink
    2014-01-05T17:44:52Z  SDO/HMI HMI_FRONT2 ...   6173.0 Invalid KeyLink
    Length = 81 rows

    You can then make the request and download the data:

    >>> res = client.fetch(response)   # doctest: +SKIP

    This returns a Results instance which can be used to watch the progress
    of the download.

    >>> res.wait(progress=True)   # doctest: +SKIP

    *Example 2*

    Query the JSOC for some AIA 171 data, and separate out the staging and the
    download steps:

    >>> import astropy.units as u
    >>> from sunpy.net import jsoc
    >>> from sunpy.net import attrs as a
    >>> client = jsoc.JSOCClient()
    >>> response = client.search(a.Time('2014/1/1T00:00:00', '2014/1/1T00:00:36'),
    ...                         a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Segment('image'),
    ...                         a.jsoc.Wavelength(171*u.AA), a.jsoc.Notify("sunpy@sunpy.org"))

    the response object holds the records that your query will return:

    >>> print(response)
            DATE         TELESCOP INSTRUME          T_OBS          WAVELNTH WAVEUNIT
    -------------------- -------- -------- ----------------------- -------- --------
    2014-01-06T15:07:12Z  SDO/AIA    AIA_3 2013-12-31T23:59:36.34Z      171 angstrom
    2014-01-06T15:07:12Z  SDO/AIA    AIA_3 2013-12-31T23:59:48.34Z      171 angstrom
    2014-01-07T15:05:10Z  SDO/AIA    AIA_3 2014-01-01T00:00:00.34Z      171 angstrom
    2014-01-07T15:05:10Z  SDO/AIA    AIA_3 2014-01-01T00:00:12.34Z      171 angstrom

    You can then make the request:

    >>> requestIDs = client.request_data(response)
    [u'JSOC_20140724_952']

    This returns a list of all the request identifiers for your query.

    You can then check the status of the request, which will print out a status
    message and return you the status code, a code of 1 means it is not ready
    to download and a code of 0 means the request is staged and ready. A code
    of 6 means an error, which is commonly that the request has not had time to
    get into the queue.

    >>> status = client.check_request(requestIDs)
    Request JSOC_20140724_955 was submitted 10 seconds ago, it is not ready to download.

    Once the status code is 0 you can download the data using the `get_request`
    method:

    >>> res = client.get_request(requestIDs)

    This returns a Results instance which can be used to watch the progress
    of the download.

    >>> res.wait(progress=True)   # doctest: +SKIP
    """

    def initialise(self, ser):

        c = drms.Client()
        pkeys = c.pkeys(ser)
        for pkey in pkeys:
            genClass = type(pkey, (_VSOSimpleAttr,), {})
            setattr(attrs, genClass.__name__, genClass)

    def search(self, *query, **kwargs):
        """
        Build a JSOC query and submit it to JSOC for processing.

        Takes a variable number of :mod:`sunpy.net.jsoc.attrs` as parameters,
        which are chained together using the AND (`&`) operator.

        Complex queries to be easily formed using logical operators such as
        `&` and `|`, in the same way as the VSO client.

        Examples
        --------
        Request all AIA 304 image data between 2010-01-01T00:00 and
        2010-01-01T01:00 in rice compressed form.

        >>> import astropy.units as u
        >>> from sunpy.net import jsoc
        >>> from sunpy.net import attrs as a
        >>> client = jsoc.JSOCClient()
        >>> response = client.search(a.Time('2010-01-01T00:00:00', '2010-01-01T01:00:00'),
        ...                         a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Wavelength(304*u.AA),
        ...                         a.jsoc.Compression('rice'), a.jsoc.Segment('image'))

        Returns
        -------
        results : JSOCResults object
            A collection of records that the query returns.
        """

        return_results = JSOCResponse()
        query = and_(*query)
        blocks = []
        for block in walker.create(query):
            iargs = kwargs.copy()
            iargs.update(block)
            blocks.append(iargs)

            return_results.append(self._lookup_records(iargs))

        return_results.query_args = blocks

        return return_results

    @deprecated('0.8', alternative='JSOCClient.search')
    def query(self, *query, **kwargs):
        """
        See `~sunpy.net.jsoc.jsoc.JSOCClient.search`
        """
        return self.search(*query, **kwargs)

    def request_data(self, jsoc_response, **kwargs):
        """
        Request that JSOC stages the data for download.

        Parameters
        ----------
        jsoc_response : JSOCResponse object
            The results of a query

        Returns
        -------
        requests : ExportRequest Object or
                   a list of ExportRequest objects
            
            Request Id can be accessed by requests.id
            Request status can be accessed by requests.status 

        """

        # Do a multi-request for each query block
        requests = []
        for block in jsoc_response.query_args:

            ds = self._make_recordset(**block)
            cd = drms.Client(email=block.get('notify', ''))
            protocol = block.get('protocol', 'fits')

            if protocol != 'fits' and protocol != 'as-is':
                error_message = "Protocols other than fits and as-is are "\
                                "are not supported."
                raise TypeError(error_message)

            method = 'url' if protocol == 'fits' else 'url-quick'
            r = cd.export(ds, method=method, protocol=protocol)
            r.wait()

            requests.append(r)

        if len(responses) == 1:
            return responses[0]
        return responses

    def check_request(self, responses):
        """
        Check the status of a request and print out a message about it

        Parameters
        ----------
        requestIDs : list or string
            A list of requestIDs to check

        Returns
        -------
        status : list
            A list of status' that were returned by JSOC
        """
        # Convert IDs to a list if not already
        if not isiterable(responses) or isinstance(responses, drms.ExportRequest):
            responses = [responses]

        allstatus = []
        for response in responses:
            status = response.status

            if status == 0:  # Data ready to download
                print("Request {0} was exported at {1} and is ready to "
                      "download.".format(response.id,
                                         response._d['exptime']))
            elif status == 1:
                print_message = "Request {0} was submitted {1} seconds ago, "\
                                "it is not ready to download."
                print(print_message.format(response.id,
                                           response._d['wait']))
            else:
                print_message = "Request returned status: {0} with error: {1}"
                json_status = response.status
                json_error = response._d['error']
                print(print_message.format(json_status, json_error))

            allstatus.append(status)

        if len(allstatus) == 1:
            return allstatus[0]
        return allstatus

    def fetch(self, jsoc_response, path=None, overwrite=False, progress=True,
            max_conn=5, downloader=None, sleep=10):
        """
        Make the request for the data in jsoc_response and wait for it to be
        staged and then download the data.

        Parameters
        ----------
        jsoc_response : JSOCResponse object
            A response object

        path : string
            Path to save data to, defaults to SunPy download dir

        overwrite : bool
            Replace files with the same name if True

        progress : bool
            Print progress info to terminal

        max_conns : int
            Maximum number of download connections.

        downloader: `sunpy.download.Downloader` instance
            A Custom downloader to use

        sleep : int
            The number of seconds to wait between calls to JSOC to check the status
            of the request.

        Returns
        -------
        results : a :class:`sunpy.net.vso.Results` instance
            A Results object
        """

        # Make staging request to JSOC
        responses = self.request_data(jsoc_response)
        # Add them to the response for good measure
        jsoc_response.requestIDs = [r.id for r in responses]
        time.sleep(sleep/2.)

        r = Results(lambda x: None, done=lambda maps: [v['path'] for v in maps.values()])

        for response in responses:

            if progress:
                self.check_request(response)

            if response.status == 0:
                res = response

                r = self.get_request(res, path=path, overwrite=overwrite,
                                     progress=progress, results=r)
            else:
                time.sleep(sleep)

        return r

    @deprecated('0.8', alternative='JSOCClient.fetch')
    def get(self, jsoc_response, path=None, overwrite=False, progress=True,
            max_conn=5, downloader=None, sleep=10):
        """
        See `~sunpy.net.jsoc.jsoc.JSOCClient.fetch`
        """
        return self.fetch(jsoc_response, path=path, overwrite=overwrite, progress=progress,
            max_conn=max_conn, downloader=downloader, sleep=sleep)


    def get_request(self, responses, path=None, overwrite=False, progress=True,
                    max_conn=5, downloader=None, results=None):
        """
        Query JSOC to see if request_id is ready for download.

        If the request is ready for download, download it.

        Parameters
        ----------
        requestIDs : list or string
            One or many requestID strings

        path : string
            Path to save data to, defaults to SunPy download dir

        overwrite : bool
            Replace files with the same name if True

        progress : bool
            Print progress info to terminal

        max_conns : int
            Maximum number of download connections.

        downloader : `sunpy.download.Downloader` instance
            A Custom downloader to use

        results: Results instance
            A Results manager to use.

        Returns
        -------
        res: Results
            A Results instance or None if no URLs to download
        """

        # Convert Responses to a list if not already

        if not isiterable(responses) or isinstance(responses, drms.ExportRequest):
            responses = [responses]

        if path is None:
            path = config.get('downloads', 'download_dir')
        path = os.path.expanduser(path)

        if downloader is None:
            downloader = Downloader(max_conn=max_conn, max_total=max_conn)

        # A Results object tracks the number of downloads requested and the
        # number that have been completed.
        if results is None:
            results = Results(lambda _: downloader.stop())
        urls = []
        for response in responses:

            if response.status == 0:
                for index, data in response.data.iterrows():

                    is_file = os.path.isfile(os.path.join(path, data['filename']))
                    if overwrite or not is_file:
                        url_dir = response.request_url + '/'
                        urls.append(urllib.parse.urljoin(url_dir, data['filename']))

                    else:
                        print_message = "Skipping download of file {} as it " \
                                        "has already been downloaded"
                        print(print_message.format(data['filename']))
                        # Add the file on disk to the output
                        results.map_.update({data['filename']:
                                            {'path': os.path.join(path, data['filename'])}})

                if progress:
                    print_message = "{0} URLs found for download. Totalling {1}MB"
                    print(print_message.format(len(urls), response._d['size']))

            else:
                if progress:
                    self.check_request(response)

        if urls:
            for url in urls:
                downloader.download(url, callback=results.require([url]),
                                    errback=lambda x: print(x), path=path)

        else:
            # Make Results think it has finished.
            results.require([])
            results.poke()

        return results

    def _process_time(self, time):
        """
        Take a UTC time string or datetime instance and generate a astropy.time
        object in TAI frame. Alternatively convert a astropy time object to TAI

        Parameters
        ----------
        time: six.string_types or datetime or astropy.time
            Input time

        Returns
        -------
        datetime, in TAI
        """

        if isinstance(time, six.string_types):
            time = parse_time(time)

        time = astropy.time.Time(time, scale='utc')
        time = time.tai  # Change the scale to TAI

        return time.datetime

    def _make_recordset(self, start_time, end_time, series, wavelength='',
                        segment='', primekeys={}, **kwargs):
        """
        Take the query arguments and build a record string.
        """

        # Extract and format segment
        if segment:
            if isinstance(segment, list):
                segment = str(segment)[1:-1].replace(' ', '').replace("'", '')
            segment = '{{{segment}}}'.format(segment=segment)

        # Extract and format sample
        sample = kwargs.get('sample', '')
        if sample:
            sample = '@{}s'.format(sample)

        # Extract and format primekeys
        pkstr = ''
        c = drms.Client()
        pkeys = c.pkeys(series)

        for pkey in pkeys:
            if pkey in ['T_OBS', 'T_REC', 'T_START']:
                pkstr += '[{start}-{end}{sample}]'.format(
                    start=start_time.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                    end=end_time.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                    sample=sample)

            elif pkey == 'WAVELNTH' and wavelength is not '':
                if isinstance(wavelength, list):
                    wavelength = [int(np.ceil(wave.to(u.AA).value)) for wave in wavelength]
                    wavelength = str(wavelength)
                else:
                    wavelength = '[{0}]'.format(int(np.ceil(wavelength.to(u.AA).value)))
                pkstr += wavelength

            elif len(primekeys) > 0:
                pkstr += '[{0}]'.format(primekeys.pop(pkey, ''))

            else:
                break

        dataset = '{series}{primekeys}{segment}'.format(series=series,
                                                        primekeys=pkstr,
                                                        segment=segment)

        return dataset

    def _lookup_records(self, iargs):
        """
        Do a LookData request to JSOC to workout what results the query returns
        """

        keywords_default = ['DATE', 'TELESCOP', 'INSTRUME', 'T_OBS', 'WAVELNTH']
        isMeta = iargs.get('meta', False)
        c = drms.Client()

        if isMeta:
            keywords = '***ALL***'
        else:
            keywords = iargs.get('keys', keywords_default)

        if not all([k in iargs for k in ('start_time', 'end_time', 'series')]):
            error_message = "Both Time and Series must be specified for a "\
                            "JSOC Query"
            raise ValueError(error_message)

        if not isinstance(keywords, list) and not isinstance(keywords, six.string_types):
            error_message = "Keywords can only be passed as a list or"\
                            " comma-separated strings."
            raise ValueError(error_message)

        segments = iargs.get('segment', '')
        if segments:
            if not isinstance(segments, list) and not isinstance(segments, six.string_types):
                error_message = "Segments can only be passed as a list or"\
                                " comma-separated strings."
                raise ValueError(error_message)

        pkeys = c.pkeys(iargs['series'])
        pkeys_passed = iargs.get('primekeys', None)
        if pkeys_passed is not None:
            if not set(list(pkeys_passed.keys())) < set(pkeys):
                error_message = "Unexpected PrimeKeys were passed. The series {series} "\
                                "supports the following PrimeKeys {pkeys}"
                raise TypeError(error_message.format(series=iargs['series'], pkeys=pkeys))

        wavelength = iargs.get('wavelength', '')
        if wavelength:
            if 'WAVELNTH' not in pkeys:
                error_message = "The series {series} does not support wavelength attribute."\
                                " Following primekeys are supported {pkeys}"
                raise TypeError(error_message.format(series=iargs['series'], pkeys=pkeys))

        si = c.info(iargs['series'])
        segs = list(si.segments.index.values)
        segs_passed = iargs.get('segment', None)
        if segs_passed is not None:
            if not set(segs_passed) < set(segs):
                error_message = "Unexpected Segments were passed. The series {series} "\
                                "contains the following Segments {segs}"
                raise TypeError(error_message.format(series=iargs['series'], segs=segs))

        iargs['start_time'] = self._process_time(iargs['start_time'])
        iargs['end_time'] = self._process_time(iargs['end_time'])

        postthis = {'ds': self._make_recordset(**iargs),
                    'key': str(keywords)[1:-1].replace(' ', '').replace("'", ''),
                    'seg': '**NONE**',
                    'link': '**NONE**',
                    }

        r = c.query(postthis['ds'], key=postthis['key'], rec_index=isMeta)

        if isMeta:
            return r

        if r is None or r.empty:
            return astropy.table.Table()
        else:
            return astropy.table.Table.from_pandas(r)

    @classmethod
    def _can_handle_query(cls, *query):
        chkattr = ['Series', 'Protocol', 'Notify', 'Wavelength', 'Time',
                   'Segment', 'Keys']

        return all([x.__class__.__name__ in chkattr for x in query])
