# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import

import os
import time
import warnings

import requests
import numpy as np
import pandas as pd
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
from sunpy.net.jsoc.attrs import walker, Keys, PrimeKey
from sunpy.net.vso.attrs import _VSOSimpleAttr
from sunpy.extern.six.moves import urllib
from sunpy.extern import six
from sunpy.util.metadata import MetaDict
from sunpy.util import deprecated

__all__ = ['JSOCClient', 'JSOCResponse']


PKEY_LIST_TIME = ['T_START', 'T_REC', 'T_OBS', 'MidTime', 'OBS_DATE',
                  'obsdate', 'DATE_OBS', 'starttime', 'stoptime', 'UTC_StartTime']


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
        data. See `this <http://jsoc.stanford.edu/ajax/register_email.html>`

    Notes
    -----
    This Client mocks input to this `site <http://jsoc.stanford.edu/ajax/exportdata.html>`
    Therefore that is a good resource if things are mis-behaving.
    The full list of 'series' is available through this `site <http://jsoc.stanford.edu>`

    You can build more complex queries by specifying parameters to POST to JSOC via keyword
    arguments. You can generate these kwargs using the Export Data page at JSOC.

    JSOC now requires a validated email address, you can pass in your validated email address
    using the `~sunpy.net.jsoc.attrs.Notify` attribute. You have to register your email address
    with JSOC beforehand `here <http://jsoc.stanford.edu/ajax/register_email.html>`.

    The backend of SunPy's JSOC Client uses `drms package <https://github.com/kbg/drms>`
    The tutorials can be found `here <https://drms.readthedocs.io/en/stable/tutorial.html>`
    This can be used to build complex queries, by directly inputting the query string.

    Examples
    --------

    *Example 1*

    Query JSOC for some HMI data at 45 second cadence:

    >>> from sunpy.net import jsoc
    >>> from sunpy.net import attrs as a
    >>> client = jsoc.JSOCClient()
    >>> response = client.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                          a.jsoc.Series('hmi.m_45s'), a.jsoc.Notify("sunpy@sunpy.org"))

    the response object holds the records that your query will return:

    >>> print(response)   # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE

    <Table length=81>
            DATE         TELESCOP  INSTRUME           T_OBS          WAVELNTH
           str20           str7     str10             str23          float64
    -------------------- -------- ---------- ----------------------- --------
    2014-01-05T17:46:02Z  SDO/HMI HMI_FRONT2 2014.01.01_00:00:37_TAI   6173.0
    2014-01-05T17:47:11Z  SDO/HMI HMI_FRONT2 2014.01.01_00:01:22_TAI   6173.0
    2014-01-05T17:48:18Z  SDO/HMI HMI_FRONT2 2014.01.01_00:02:07_TAI   6173.0
    2014-01-05T17:49:26Z  SDO/HMI HMI_FRONT2 2014.01.01_00:02:52_TAI   6173.0
    2014-01-05T17:50:34Z  SDO/HMI HMI_FRONT2 2014.01.01_00:03:37_TAI   6173.0
    2014-01-05T17:51:43Z  SDO/HMI HMI_FRONT2 2014.01.01_00:04:22_TAI   6173.0
                     ...      ...        ...                     ...      ...
    2014-01-05T17:40:18Z  SDO/HMI HMI_FRONT2 2014.01.01_00:56:52_TAI   6173.0
    2014-01-05T17:41:25Z  SDO/HMI HMI_FRONT2 2014.01.01_00:57:37_TAI   6173.0
    2014-01-05T17:42:33Z  SDO/HMI HMI_FRONT2 2014.01.01_00:58:22_TAI   6173.0
    2014-01-05T17:43:41Z  SDO/HMI HMI_FRONT2 2014.01.01_00:59:07_TAI   6173.0
    2014-01-05T17:44:52Z  SDO/HMI HMI_FRONT2 2014.01.01_00:59:52_TAI   6173.0
    2014-01-05T17:46:04Z  SDO/HMI HMI_FRONT2 2014.01.01_01:00:37_TAI   6173.0

    Length = 81 rows

    You can then make the request and download the data:

    >>> res = client.fetch(response)   # doctest: +SKIP

    This returns a Results instance which can be used to watch the progress
    of the download.

    Note
    ----
    A registered email address is not required if you only need to query for data,
    it is used only if you need to make an export request. For example,

    >>> client = jsoc.JSOCClient()
    >>> response = client.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                          a.jsoc.Series('hmi.m_45s'))

    The above is a successful query operation, and will return query responses as before.

    But, this response object cannot be used to make an export request and will throw an
    error if done so:

    >>> res = client.fetch(response)   # doctest: +SKIP

    ValueError: Email address is invalid or not registered

    *Example 2*

    Query the JSOC for some AIA 171 data, and separate out the staging and the
    download steps:

    >>> import astropy.units as u
    >>> from sunpy.net import jsoc
    >>> from sunpy.net import attrs as a
    >>> client = jsoc.JSOCClient()
    >>> response = client.search(a.jsoc.Time('2014/1/1T00:00:00', '2014/1/1T00:00:36'),
    ...                          a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Segment('image'),
    ...                          a.jsoc.Wavelength(171*u.AA), a.jsoc.Notify("sunpy@sunpy.org"))

    the response object holds the records that your query will return:

    >>> print(response)

    <Table length=4>
            DATE         TELESCOP INSTRUME          T_OBS          WAVELNTH
           str20           str7     str5            str23           int64
    -------------------- -------- -------- ----------------------- --------
    2014-01-07T15:05:10Z  SDO/AIA    AIA_3 2014-01-01T00:00:12.34Z      171
    2014-01-07T15:05:10Z  SDO/AIA    AIA_3 2014-01-01T00:00:24.34Z      171
    2014-01-07T15:05:10Z  SDO/AIA    AIA_3 2014-01-01T00:00:36.34Z      171
    2014-01-07T15:05:10Z  SDO/AIA    AIA_3 2014-01-01T00:00:48.34Z      171


    You can then make the request:

    >>> requests = client.request_data(response)
    <ExportRequest id="JSOC_20170713_1461", status=0>

    This returns a list of all the ExportRequest objects for your query. You can
    get the ExportRequest ID :

    >>> requests.id
    'JSOC_20170713_1461'

    >>> requests.status
    0

    You can also check the status of the request, which will print out a status
    message and return you the status code, a code of 1 means it is not ready
    to download and a code of 0 means the request is staged and ready. A code
    of 6 means an error, which is commonly that the request has not had time to
    get into the queue.

    >>> requests.status
    0

    or

    >>> status = client.check_request(requests)
    Request JSOC_20140724_955 was submitted 10 seconds ago, it is not ready to download.

    Once the status code is 0 you can download the data using the `get_request`
    method:

    >>> res = client.get_request(requests)

    This returns a Results instance which can be used to watch the progress
    of the download.

    >>> res.wait(progress=True)   # doctest: +SKIP
    """

    def search(self, *query, **kwargs):
        """
        Build a JSOC query and submit it to JSOC for processing.

        Takes a variable number of :mod:`sunpy.net.jsoc.attrs` as parameters,
        which are chained together using the AND (`&`) operator.

        Complex queries to be easily formed using logical operators such as
        `&` and `|`, in the same way as the VSO client.

        Parameters
        ----------
        query : a variable number of :mod:`sunpy.net.jsoc.attrs`
                as parameters, which are chained together using
                the AND (`&`) operator.

        Returns
        -------
        response : `~sunpy.net.jsoc.jsoc.JSOCResponse` object
            A collection of records that the query returns.

        Examples
        --------

        *Example 1*

        Request all AIA 304 image data between 2014-01-01T00:00 and
        2014-01-01T01:00.

        >>> import astropy.units as u
        >>> from sunpy.net import jsoc
        >>> from sunpy.net import attrs as a
        >>> client = jsoc.JSOCClient()
        >>> response = client.search(a.jsoc.Time('2010-01-01T00:00:00', '2010-01-01T01:00:00'),
        ...                          a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Wavelength(304*u.AA),
        ...                          a.jsoc.Segment('image'))


        *Example 2*

        Request keyword data of hmi.v_45s for certain specific keywords only.

        >>> import astropy.units as u
        >>> from sunpy.net import jsoc
        >>> from sunpy.net import attrs as a
        >>> client = jsoc.JSOCClient()
        >>> response = client.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
        ...                          a.jsoc.Series('hmi.v_45s'),
        ...                          a.jsoc.Keys('T_REC, DATAMEAN, OBS_VR'))
        >>> print(response)

        <Table length=81>
                 T_REC            DATAMEAN     OBS_VR
                 str23            float64     float64
        ----------------------- ----------- -----------
        2014.01.01_00:00:45_TAI 1906.518188 1911.202614
        2014.01.01_00:01:30_TAI 1908.876221 1913.945512
        2014.01.01_00:02:15_TAI   1911.7771 1916.667999
        2014.01.01_00:03:00_TAI 1913.422485 1919.369924
        2014.01.01_00:03:45_TAI 1916.500488 1922.050862
                            ...         ...         ...
        2014.01.01_00:57:45_TAI 2054.584473 2058.971861
        2014.01.01_00:58:30_TAI 2056.094238 2060.075964
        2014.01.01_00:59:15_TAI 2056.366699 2061.157734
        2014.01.01_01:00:00_TAI 2057.013428 2062.217153
        2014.01.01_01:00:45_TAI 2059.014893 2063.254285

        *Example 3*

        Request data of aia.lev1_euv_12s on the basis of primekeys other than T_REC.

        >>> import astropy.units as u
        >>> from sunpy.net import jsoc
        >>> from sunpy.net import attrs as a
        >>> client = jsoc.JSOCClient()
        >>> response = client.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T02:00:00'),
        ...                          a.jsoc.Series('aia.lev1_euv_12s'),
        ...                          a.jsoc.PrimeKey('WAVELNTH','171'))

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

    def search_metadata(self, *query, **kwargs):
        """
        Get the metadata of all the files obtained in a search query.
        Builds a jsoc query, similar to query() method, and takes similar inputs.

        Complex queries to be easily formed using logical operators such as
        `&` and `|`, in the same way as the query() function.

        Parameters
        ----------
        query : a variable number of :mod:`sunpy.net.jsoc.attrs`
                as parameters, which are chained together using
                the AND (`&`) operator.

        Returns
        -------
        res : `~pandas.DataFrame` object
            A collection of metadata of all the files.

        Example
        -------

        Request metadata or all all AIA 304 image data between 2014-01-01T00:00 and
        2014-01-01T01:00.

        Since, the function only performs a lookdata, and does not make a proper export
        request, attributes like Segment need not be passed.

        >>> import astropy.units as u
        >>> from sunpy.net import jsoc
        >>> from sunpy.net import attrs as a
        >>> client = jsoc.JSOCClient()
        >>> metadata = client.search_metadata(
                                    a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T00:02:00'),
        ...                         a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Wavelength(304*u.AA))
        >>> print(metadata)

                                                          T_REC                  ...         T_REC_epoch
        aia.lev1_euv_12s[2014-01-01T00:00:01Z][304]  2014-01-01T00:00:01Z        ...     1993.01.01_00:00:04_TAI
        aia.lev1_euv_12s[2014-01-01T00:00:13Z][304]  2014-01-01T00:00:13Z        ...     1993.01.01_00:00:04_TAI
        aia.lev1_euv_12s[2014-01-01T00:00:25Z][304]  2014-01-01T00:00:25Z        ...     1993.01.01_00:00:04_TAI
        aia.lev1_euv_12s[2014-01-01T00:00:37Z][304]  2014-01-01T00:00:37Z        ...     1993.01.01_00:00:04_TAI
        aia.lev1_euv_12s[2014-01-01T00:00:49Z][304]  2014-01-01T00:00:49Z        ...     1993.01.01_00:00:04_TAI
        aia.lev1_euv_12s[2014-01-01T00:01:01Z][304]  2014-01-01T00:01:01Z        ...     1993.01.01_00:00:04_TAI
        aia.lev1_euv_12s[2014-01-01T00:01:13Z][304]  2014-01-01T00:01:13Z        ...     1993.01.01_00:00:04_TAI
        aia.lev1_euv_12s[2014-01-01T00:01:25Z][304]  2014-01-01T00:01:25Z        ...     1993.01.01_00:00:04_TAI
        aia.lev1_euv_12s[2014-01-01T00:01:37Z][304]  2014-01-01T00:01:37Z        ...     1993.01.01_00:00:04_TAI
        aia.lev1_euv_12s[2014-01-01T00:01:49Z][304]  2014-01-01T00:01:49Z        ...     1993.01.01_00:00:04_TAI
        aia.lev1_euv_12s[2014-01-01T00:02:01Z][304]  2014-01-01T00:02:01Z        ...     1993.01.01_00:00:04_TAI


        [11 rows x 176 columns]

        """
        query = and_(*query)
        blocks = []
        res = pd.DataFrame()
        for block in walker.create(query):
            iargs = kwargs.copy()
            iargs.update(block)
            iargs.update({'meta': True})
            blocks.append(iargs)
            res = res.append(self._lookup_records(iargs))

        return res

    def request_data(self, jsoc_response, **kwargs):
        """
        Request that JSOC stages the data for download.

        Parameters
        ----------
        jsoc_response : `~sunpy.net.jsoc.jsoc.JSOCResponse` object
            The results of a query

        Returns
        -------
        requests : `~drms.ExportRequest` Object or
                   a list of  `~drms.ExportRequest` objects

                   Request Id can be accessed by requests.id
                   Request status can be accessed by requests.status

        """

        requests = []
        for block in jsoc_response.query_args:

            ds = self._make_recordset(**block)
            cd = drms.Client(email=block.get('notify', ''))
            protocol = block.get('protocol', 'fits')

            if protocol != 'fits' and protocol != 'as-is':
                error_message = "Protocols other than fits and as-is are "\
                                "are not supported."
                raise TypeError(error_message)

            method = 'url' if protocol == 'fits' else 'url_quick'
            r = cd.export(ds, method=method, protocol=protocol)
            r.wait()

            requests.append(r)

        if len(requests) == 1:
            return requests[0]
        return requests

    def check_request(self, requests):
        """
        Check the status of a request and print out a message about it

        Parameters
        ----------
        requests :  `~drms.ExportRequest` object or
                   a list of  `~drms.ExportRequest` objects,
                   returned by `~sunpy.net.jsoc.jsoc.JSOCClient.request_data`

        Returns
        -------
        status : int or list
            A status or list of status' that were returned by JSOC
        """
        # Convert IDs to a list if not already
        if not isiterable(requests) or isinstance(requests, drms.ExportRequest):
            requests = [requests]

        allstatus = []
        for request in requests:
            status = request.status

            if status == 0:  # Data ready to download
                print("Request {0} was exported at {1} and is ready to "
                      "download.".format(request.id,
                                         request._d['exptime']))
            elif status == 1:
                print_message = "Request {0} was submitted {1} seconds ago, "\
                                "it is not ready to download."
                print(print_message.format(request.id,
                                           request._d['wait']))
            else:
                print_message = "Request returned status: {0} with error: {1}"
                json_status = request.status
                json_error = request._d['error']
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
        jsoc_response : `~sunpy.net.jsoc.jsoc.JSOCResponse` object
            A response object

        path : string
            Path to save data to, defaults to SunPy download dir

        overwrite : bool
            Replace files with the same name if True

        progress : bool
            Print progress info to terminal

        max_conns : int
            Maximum number of download connections.

        downloader: `~sunpy.download.Downloader` instance
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
        # Make response iterable
        if not isiterable(responses):
            responses = [responses]
        # Add them to the response for good measure
        jsoc_response.requestIDs = [r.id for r in responses]
        time.sleep(sleep/2.)

        r = Results(lambda x: None, done=lambda maps: [v['path'] for v in maps.values()])

        for response in responses:

            if progress:
                self.check_request(response)
                r = self.get_request(response, path=path, overwrite=overwrite,
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

    def get_request(self, requestIDs, path=None, overwrite=False, progress=True,
                    max_conn=5, downloader=None, results=None):
        """
        Query JSOC to see if request_id is ready for download.

        If the request is ready for download, download it.

        Parameters
        ----------
        requests : `~drms.ExportRequest` object or
                   a list of `~drms.ExportRequest` objects,
                   returned by `~sunpy.net.jsoc.jsoc.JSOCClient.request_data`

        path : string
            Path to save data to, defaults to SunPy download dir

        overwrite : bool
            Replace files with the same name if True

        progress : bool
            Print progress info to terminal

        max_conns : int
            Maximum number of download connections.

        downloader : `~sunpy.download.Downloader` instance
            A Custom downloader to use

        results: Results instance
            A Results manager to use.

        Returns
        -------
        res: Results
            A Results instance or None if no URLs to download
        """

        # Convert Responses to a list if not already

        if not isiterable(responses):
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

    def _make_recordset(self, series, start_time='', end_time='', wavelength='',
                        segment='', primekey={}, **kwargs):
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

        # Populate primekeys dict with Time and Wavelength values
        if start_time and end_time:
            if not any(x in PKEY_LIST_TIME for x in primekey):
                timestr = '{start}-{end}{sample}'.format(
                        start=start_time.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                        end=end_time.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                        sample=sample)
            else:
                error_message = "Time attribute has been passed both as a Time()"\
                                " and PrimeKey(). Please provide any one of them"\
                                " or separate them by OR operator."
                raise ValueError(error_message)

        else:
            for i in PKEY_LIST_TIME:
                timestr = '{0}'.format(primekey.pop(i, ''))
                if timestr:
                    break

        if wavelength:
            if not primekey.get('WAVELNTH', ''):
                if isinstance(wavelength, list):
                    wavelength = [int(np.ceil(wave.to(u.AA).value)) for wave in wavelength]
                    wavelength = str(wavelength)
                else:
                    wavelength = '{0}'.format(int(np.ceil(wavelength.to(u.AA).value)))

            else:
                error_message = "Wavelength attribute has been passed both as a Wavelength()"\
                                " and PrimeKey(). Please provide any one of them"\
                                " or separate them by OR operator."
                raise ValueError(error_message)

        else:
            wavelength = '{0}'.format(primekey.pop('WAVELNTH', ''))

        if timestr:
            primekey['TIME'] = timestr
        if wavelength:
            primekey['WAVELNTH'] = wavelength

        # Extract and format primekeys
        pkstr = ''
        c = drms.Client()
        si = c.info(series)
        pkeys_isTime = si.keywords.loc[si.primekeys].is_time
        for pkey in pkeys_isTime.index.values:

            if len(primekey) > 0:
                if pkeys_isTime[pkey]:
                    pkstr += '[{0}]'.format(primekey.pop('TIME', ''))
                else:
                    pkstr += '[{0}]'.format(primekey.pop(pkey, ''))
            else:
                break

        if not pkstr:
            error_message = "Atleast one PrimeKey must be passed."
            raise ValueError(error_message)

        dataset = '{series}{primekey}{segment}'.format(series=series,
                                                       primekey=pkstr,
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
            keywords = '**ALL**'
        else:
            keywords = iargs.get('keys', keywords_default)

        if 'series' not in iargs:
            error_message = "Series must be specified for a JSOC Query"
            raise ValueError(error_message)

        if not isinstance(keywords, list) and not isinstance(keywords, six.string_types):
            error_message = "Keywords can only be passed as a list or "\
                            "comma-separated strings."
            raise TypeError(error_message)

        # Raise errors for PrimeKeys
        pkeys = c.pkeys(iargs['series'])
        pkeys_passed = iargs.get('primekey', None)
        if pkeys_passed is not None:
            if not set(list(pkeys_passed.keys())) <= set(pkeys):
                error_message = "Unexpected PrimeKeys were passed. The series {series} "\
                                "supports the following PrimeKeys {pkeys}"
                raise ValueError(error_message.format(series=iargs['series'], pkeys=pkeys))

        # Raise errors for wavelength
        wavelength = iargs.get('wavelength', '')
        if wavelength:
            if 'WAVELNTH' not in pkeys:
                error_message = "The series {series} does not support wavelength attribute."\
                                "The following primekeys are supported {pkeys}"
                raise TypeError(error_message.format(series=iargs['series'], pkeys=pkeys))

        # Raise errors for segments
        si = c.info(iargs['series'])
        segs = list(si.segments.index.values)
        segs_passed = iargs.get('segment', '')
        if segs_passed:
            if not isinstance(segs_passed, list):
                error_message = "Segments can only be passed as a list of strings."
                raise TypeError(error_message)

            if not set(segs_passed) <= set(segs):
                error_message = "Unexpected Segments were passed. The series {series} "\
                                "contains the following Segments {segs}"
                raise ValueError(error_message.format(series=iargs['series'], segs=segs))

        if 'start_time' in iargs:
            iargs['start_time'] = iargs['start_time'].tai.datetime
            iargs['end_time'] = iargs['end_time'].tai.datetime

        ds = self._make_recordset(**iargs)
        if isinstance(keywords, list):
            key = str(keywords)[1:-1].replace(' ', '').replace("'", '')
        else:
            key = keywords

        r = c.query(ds, key=key, rec_index=isMeta)

        if isMeta:
            return r

        if r is None or r.empty:
            return astropy.table.Table()
        else:
            return astropy.table.Table.from_pandas(r)

    @classmethod
    def _can_handle_query(cls, *query):
        chkattr = ['Series', 'Protocol', 'Notify', 'Wavelength', 'Time',
                   'Segment', 'Keys', 'PrimeKey']

        return all([x.__class__.__name__ in chkattr for x in query])
