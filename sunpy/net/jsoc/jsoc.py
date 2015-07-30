# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import

import os
import time
import urlparse
import warnings

import requests
import numpy as np
import astropy.units as u
import astropy.time
import astropy.table
from astropy.utils.misc import isiterable

from sunpy import config
from sunpy.time import parse_time, TimeRange
from sunpy.net.download import Downloader
from sunpy.net.vso.vso import Results
from sunpy.net.attr import and_
from sunpy.net.jsoc.attrs import walker

__all__ = ['JSOCClient', 'JSOCResponse']

JSOC_INFO_URL = 'http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_info'
JSOC_EXPORT_URL = 'http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_fetch'
BASE_DL_URL = 'http://jsoc.stanford.edu'


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
    >>> client = jsoc.JSOCClient()
    >>> response = client.query(jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                         jsoc.Series('hmi.m_45s'), jsoc.Notify("sunpy@sunpy.org"))

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

    >>> res = client.get(response)   # doctest: +SKIP

    This returns a Results instance which can be used to watch the progress
    of the download.

    >>> res.wait(progress=True)   # doctest: +SKIP

    *Example 2*

    Query the JSOC for some AIA 171 data, and separate out the staging and the
    download steps:

    >>> import astropy.units as u
    >>> from sunpy.net import jsoc
    >>> client = jsoc.JSOCClient()
    >>> response = client.query(jsoc.Time('2014/1/1T00:00:00', '2014/1/1T00:00:36'),
    ...                         jsoc.Series('aia.lev1_euv_12s'), jsoc.Segment('image'),
    ...                         jsoc.Wavelength(171*u.AA), jsoc.Notify("sunpy@sunpy.org"))

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

    This returns a list of all the request identifiers for your query.

    You can then check the status of the request, which will print out a status
    message and return you the status code, a code of 1 means it is not ready
    to download and a code of 0 means the request is staged and ready. A code
    of 6 means an error, which is commonly that the request has not had time to
    get into the queue.

    >>> status = client.check_request(requestIDs)

    Once the status code is 0 you can download the data using the `get_request`
    method:

    >>> res = client.get_request(requestIDs)

    This returns a Results instance which can be used to watch the progress
    of the download.

    >>> res.wait(progress=True)   # doctest: +SKIP
    """

    def query(self, *query, **kwargs):
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
        >>> client = jsoc.JSOCClient()
        >>> response = client.query(jsoc.Time('2010-01-01T00:00:00', '2010-01-01T01:00:00'),
        ...                         jsoc.Series('aia.lev1_euv_12s'), jsoc.Wavelength(304*u.AA),
        ...                         jsoc.Compression('rice'), jsoc.Segment('image'))

        Returns
        -------
        results : JSOCResults object
            A collection of records that the query returns.
        """

        return_results = JSOCResponse()
        query = and_(*query)
        for block in walker.create(query):
            iargs = kwargs.copy()
            iargs.update(block)

            return_results.append(self._lookup_records(iargs))

        return_results.query_args = iargs

        return return_results

    def request_data(self, jsoc_response, **kwargs):
        """
        Request that JSOC stages the data for download.

        Parameters
        ----------
        jsoc_response : JSOCResponse object
            The results of a query

        Returns
        -------
        requestIDs : list of strings
            List of the JSOC request identifiers

        """
        # A little (hidden) debug feature
        return_responses = kwargs.pop('return_resp', False)
        if len(kwargs):
            warn_message = "request_data got unexpected keyword arguments {0}"
            raise TypeError(warn_message.format(kwargs.keys()))

        # Do a multi-request for each query block
        responses = self._multi_request(**jsoc_response.query_args)
        for i, response in enumerate(responses):
            if response.status_code != 200:
                warn_message = "Query {0} retuned code {1}"
                warnings.warn(
                    Warning(warn_message.format(i, response.status_code)))
                responses.pop(i)
            elif response.json()['status'] != 2:
                warn_message = "Query {0} returned status {1} with error {2}"
                json_response = response.json()
                json_status = json_response['status']
                json_error = json_response['error']
                warnings.warn(Warning(warn_message.format(i, json_status,
                                                          json_error)))
                responses.pop(i)

        # Extract the IDs from the JSON
        requestIDs = [response.json()['requestid'] for response in responses]

        if return_responses:
            return responses

        return requestIDs

    def check_request(self, requestIDs):
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
        if not isiterable(requestIDs) or isinstance(requestIDs, basestring):
            requestIDs = [requestIDs]

        allstatus = []
        for request_id in requestIDs:
            u = self._request_status(request_id)
            status = int(u.json()['status'])

            if status == 0:  # Data ready to download
                print("Request {0} was exported at {1} and is ready to "\
                      "download.".format(u.json()['requestid'],
                                           u.json()['exptime']))
            elif status == 1:
                print_message = "Request {0} was submitted {1} seconds ago, " \
                                "it is not ready to download."
                print(print_message.format(u.json()['requestid'],
                                           u.json()['wait']))
            else:
                print_message = "Request returned status: {0} with error: {1}"
                json_status = u.json()['status']
                json_error = u.json()['error']
                print(print_message.format(json_status, json_error))

            allstatus.append(status)

        return allstatus

    def get(self, jsoc_response, path=None, overwrite=False, progress=True,
            max_conn=5, downloader=None,sleep=10):
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
        requestIDs = self.request_data(jsoc_response)
        # Add them to the response for good measure
        jsoc_response.requestIDs = requestIDs
        time.sleep(sleep/2.)

        r = None   # Needed when all requests ends in the else branch
        while requestIDs:
            for i, request_id in enumerate(requestIDs):
                u = self._request_status(request_id)

                if progress:
                    self.check_request(request_id)

                if u.status_code == 200 and u.json()['status'] == '0':
                    rID = requestIDs.pop(i)
                    r = self.get_request(rID, path=path, overwrite=overwrite,
                                         progress=progress,downloader=downloader)

                else:
                    time.sleep(sleep)

        return r

    def get_request(self, requestIDs, path=None, overwrite=False, progress=True,
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

        # Convert IDs to a list if not already

        if not isiterable(requestIDs) or isinstance(requestIDs, basestring):
            requestIDs = [requestIDs]

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
        for request_id in requestIDs:
            u = self._request_status(request_id)

            if u.status_code == 200 and u.json()['status'] == '0':
                for ar in u.json()['data']:
                    is_file = os.path.isfile(os.path.join(path, ar['filename']))
                    if overwrite or not is_file:
                        url_dir = BASE_DL_URL + u.json()['dir'] + '/'
                        urls.append(urlparse.urljoin(url_dir, ar['filename']))

                    else:
                        print_message = "Skipping download of file {} as it " \
                                        "has already been downloaded"
                        print(print_message.format(ar['filename']))
                        # Add the file on disk to the output
                        results.map_.update({ar['filename']:{'path':os.path.join(path, ar['filename'])}})

                if progress:
                    print_message = "{0} URLs found for download. Totalling {1}MB"
                    print(print_message.format(len(urls), u.json()['size']))

            else:
                if progress:
                    self.check_request(request_id)

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
        time: basestring or datetime or astropy.time
            Input time

        Returns
        -------
        datetime, in TAI
        """
        # Convert from any input (in UTC) to TAI
        if isinstance(time, basestring):
            time = parse_time(time)

        time = astropy.time.Time(time, scale='utc')
        time = time.tai  # Change the scale to TAI

        return time.datetime

    def _make_recordset(self, start_time, end_time, series, wavelength='',
                        segment='', **kwargs):
        # Build the dataset string
        # Extract and format Wavelength
        if wavelength:
            if not series.startswith('aia'):
                raise TypeError("This series does not support the wavelength attribute.")
            else:
                if isinstance(wavelength, list):
                    wavelength = [int(np.ceil(wave.to(u.AA).value)) for wave in wavelength]
                    wavelength = str(wavelength)
                else:
                    wavelength = '[{0}]'.format(int(np.ceil(wavelength.to(u.AA).value)))

        # Extract and format segment
        if segment != '':
            segment = '{{{segment}}}'.format(segment=segment)

        sample = kwargs.get('sample', '')
        if sample:
            sample = '@{}s'.format(sample)

        dataset = '{series}[{start}-{end}{sample}]{wavelength}{segment}'.format(
                   series=series, start=start_time.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                   end=end_time.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                   sample=sample,
                   wavelength=wavelength, segment=segment)

        return dataset

    def _make_query_payload(self, start_time, end_time, series, notify=None,
                            protocol='FITS', compression='rice', **kwargs):
        """
        Build the POST payload for the query parameters
        """

        if protocol.upper() == 'FITS' and compression and compression.lower() == 'rice':
            jprotocol = 'FITS,compress Rice'
        elif protocol.upper() == 'FITS':
            jprotocol = 'FITS, **NONE**'
        else:
            jprotocol = protocol

        if not notify:
            raise ValueError("JSOC queries now require a valid email address "
                             "before they will be accepted by the server")

        dataset = self._make_recordset(start_time, end_time, series, **kwargs)
        kwargs.pop('wavelength', None)
        kwargs.pop('sample',None)

        # Build full POST payload
        payload = {'ds': dataset,
                   'format': 'json',
                   'method': 'url',
                   'notify': notify,
                   'op': 'exp_request',
                   'process': 'n=0|no_op',
                   'protocol': jprotocol,
                   'requestor': 'none',
                   'filenamefmt': '{0}.{{T_REC:A}}.{{CAMERA}}.{{segment}}'.format(series)}

        payload.update(kwargs)
        return payload

    def _send_jsoc_request(self, start_time, end_time, series, notify=None,
                           protocol='FITS', compression='rice', **kwargs):
        """
        Request that JSOC stages data for download

        This routine puts in a POST request to JSOC
        """

        payload = self._make_query_payload(start_time, end_time, series,
                                           notify=notify, protocol=protocol,
                                           compression=compression, **kwargs)

        r = requests.post(JSOC_EXPORT_URL, data=payload)

        if r.status_code != 200:
            exception_message = "JSOC POST Request returned code {0}"
            raise Exception(exception_message.format(r.status_code))

        return r, r.json()

    def _lookup_records(self, iargs):
        """
        Do a LookData request to JSOC to workout what results the query returns
        """
        keywords = ['DATE', 'TELESCOP', 'INSTRUME', 'T_OBS', 'WAVELNTH',
                    'WAVEUNIT']

        if not all([k in iargs for k in ('start_time', 'end_time', 'series')]):
            error_message = "Both Time and Series must be specified for a "\
                            "JSOC Query"
            raise ValueError(error_message)

        postthis = {'ds': self._make_recordset(**iargs),
                    'op': 'rs_list',
                    'key': str(keywords)[1:-1].replace(' ', '').replace("'", ''),
                    'seg': '**NONE**',
                    'link': '**NONE**'}

        r = requests.get(JSOC_INFO_URL, params=postthis)

        result = r.json()

        out_table = {}
        if 'keywords' in result:
            for col in result['keywords']:
                out_table.update({col['name']:col['values']})

            # sort the table before returning
            return astropy.table.Table(out_table)[keywords]

        else:
            return astropy.table.Table()

    def _multi_request(self, **kwargs):
        """
        Make a series of requests to avoid the 100GB limit
        """
        start_time = kwargs.pop('start_time', None)
        end_time = kwargs.pop('end_time', None)
        series = kwargs.pop('series', None)
        if any(x is None for x in (start_time, end_time, series)):
            return []
        start_time = self._process_time(start_time)
        end_time = self._process_time(end_time)
        tr = TimeRange(start_time, end_time)
        returns = []
        response, json_response = self._send_jsoc_request(start_time, end_time,\
                                                          series, **kwargs)

        # We skip these lines because a massive request is not a practical test.
        error_response = 'Request exceeds max byte limit of 100000MB'
        if (json_response['status'] == 3 and
            json_response['error'] == error_response):  # pragma: no cover
            returns.append(self._multi_request(tr.start(), tr.center(), series, **kwargs)[0])  # pragma: no cover
            returns.append(self._multi_request(tr.center(), tr.end(), series, **kwargs)[0])  # pragma: no cover
        else:
            returns.append(response)

        return returns

    def _request_status(self, request_id):
        """
        GET the status of a request ID
        """
        payload = {'op':'exp_status', 'requestid':request_id}
        u = requests.get(JSOC_EXPORT_URL, params=payload)

        return u

