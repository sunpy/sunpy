# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import

import os
import time
import warnings
from functools import partial
from collections import Sequence

import numpy as np
import pandas as pd
import astropy.units as u
import astropy.time
import astropy.table
from astropy.utils.misc import isiterable
import drms

from sunpy import config
from sunpy.net.download import Downloader, Results
from sunpy.net.attr import and_
from sunpy.net.jsoc.attrs import walker
from sunpy.extern.six.moves import urllib
from sunpy.extern import six
from sunpy.util import deprecated

__all__ = ['JSOCClient', 'JSOCResponse']


PKEY_LIST_TIME = {'T_START', 'T_REC', 'T_OBS', 'MidTime', 'OBS_DATE',
                  'obsdate', 'DATE_OBS', 'starttime', 'stoptime', 'UTC_StartTime'}


def simple_path(path, sock, url):
    return path


class NotExportedError(Exception):
    pass


class JSOCResponse(Sequence):
    def __init__(self, table=None):
        """
        table : `astropy.table.Table`
        """

        self.table = table
        self.query_args = None
        self.requests = None

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

    def __getitem__(self, item):
        return type(self)(self.table[item])

    def __iter__(self):
        return (t for t in [self])

    def append(self, table):
        if self.table is None:
            self.table = table
        else:
            self.table = astropy.table.vstack([self.table, table])

    def response_block_properties(self):
        """
        Returns a set of class attributes on all the response blocks.
        """
        warnings.warn("The JSOC client does not support response block properties", UserWarning)
        return set()


class JSOCClient(object):
    """
    This is a Client to the JSOC Data Export service.

    It exposes a similar API to the VSO client, although the underlying model
    is more complex. The JSOC stages data before you can download it, so a JSOC
    query is a three stage process. First you query the JSOC for records,
    a table of these records is returned. Then you can request these records to
    be staged for download and then you can download them.
    The last two stages of this process are bundled together into the `fetch()`
    method, but they can be separated if you are performing a large or complex
    query.

    .. warning::

        JSOC requires you to register your email address before requesting
        data. `See this on how to register <http://jsoc.stanford.edu/ajax/register_email.html>`__.

    Notes
    -----
    The full list of ``Series`` is available through this `site <http://jsoc.stanford.edu>`_.

    JSOC requires a validated email address, you can pass in your validated email address
    using the `~sunpy.net.jsoc.attrs.Notify` attribute. You have to register your email address
    with JSOC beforehand `here <http://jsoc.stanford.edu/ajax/register_email.html>`_.

    The backend of SunPy's JSOC Client uses `drms package <https://github.com/sunpy/drms>`_.
    The tutorials can be `found here <http://docs.sunpy.org/projects/en/stable/tutorial.html>`_.
    This can be used to build complex queries, by directly inputting the query string.

    Examples
    --------

    *Example 1*

    Query JSOC for some HMI data at 45 second cadence::

        >>> from sunpy.net import jsoc
        >>> from sunpy.net import attrs as a
        >>> client = jsoc.JSOCClient()
        >>> response = client.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T00:10:00'),
        ...                          a.jsoc.Series('hmi.m_45s'), a.jsoc.Notify("sunpy@sunpy.org"))  # doctest: +REMOTE_DATA

        The response object holds the records that your query will return:

        >>> print(response)   # doctest: +ELLIPSIS  +REMOTE_DATA
                T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
        ----------------------- -------- ---------- -------- -------
        2014.01.01_00:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:01:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:02:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:03:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:03:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:04:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:05:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:06:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:06:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:07:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:08:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:09:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:09:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
        2014.01.01_00:10:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145

    You can then make the request and download the data::

        >>> res = client.fetch(response)   # doctest: +SKIP

    This returns a Results instance which can be used to watch the progress
    of the download.

    Note
    ----
    A registered email address is not required if you only need to query for data,
    it is used only if you need to make an export request. For example,::

        >>> client = jsoc.JSOCClient()  # doctest: +REMOTE_DATA
        >>> response = client.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T00:10:00'),
        ...                          a.jsoc.Series('hmi.m_45s'))  # doctest: +REMOTE_DATA

    The above is a successful query operation, and will return query responses as before.

    But, this response object cannot be used to make an export request and will throw an
    error if done so::

        >>> res = client.fetch(response)   # doctest: +SKIP

        ValueError: Email address is invalid or not registered


    *Example 2*

    Query the JSOC for some AIA 171 data, and separate out the staging and the
    download steps::

        >>> import astropy.units as u
        >>> from sunpy.net import jsoc
        >>> from sunpy.net import attrs as a
        >>> client = jsoc.JSOCClient()  # doctest: +REMOTE_DATA
        >>> response = client.search(a.jsoc.Time('2014/1/1T00:00:00', '2014/1/1T00:00:36'),
        ...                          a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Segment('image'),
        ...                          a.jsoc.Wavelength(171*u.AA), a.jsoc.Notify("sunpy@sunpy.org"))  # doctest: +REMOTE_DATA

        The response object holds the records that your query will return:

        >>> print(response)  # doctest: +REMOTE_DATA
               T_REC         TELESCOP INSTRUME WAVELNTH CAR_ROT
        -------------------- -------- -------- -------- -------
        2014-01-01T00:00:01Z  SDO/AIA    AIA_3      171    2145
        2014-01-01T00:00:13Z  SDO/AIA    AIA_3      171    2145
        2014-01-01T00:00:25Z  SDO/AIA    AIA_3      171    2145
        2014-01-01T00:00:37Z  SDO/AIA    AIA_3      171    2145

    You can then make the request::

        >>> requests = client.request_data(response)  # doctest: +SKIP

    This returns a list of all the ExportRequest objects for your query. You can
    get the ExportRequest ID ::

        >>> requests.id  # doctest: +SKIP
        'JSOC_20171205_372'

    You can also check the status of the request, which will print out a status
    message and return you the status code, a code of 1 means it is not ready
    to download and a code of 0 means the request is staged and ready. A code
    of 6 means an error, which is commonly that the request has not had time to
    get into the queue::

        >>> requests.status  # doctest: +SKIP
        0

    Once the status code is 0 you can download the data using the `get_request`
    method::

        >>> res = client.get_request(requests)  # doctest: +SKIP

    This returns a Results instance which can be used to watch the progress
    of the download::

        >>> res.wait(progress=True)   # doctest: +SKIP

    """

    def search(self, *query, **kwargs):
        """
        Build a JSOC query and submit it to JSOC for processing.

        Takes a variable number of `~sunpy.net.jsoc.attrs` as parameters,
        which are chained together using the AND (`&`) operator.

        Complex queries to be easily formed using logical operators such as
        `&` and `|`, in the same way as the VSO client.

        Parameters
        ----------
        query : a variable number of `~sunpy.net.jsoc.attrs`
                as parameters, which are chained together using
                the ``AND`` (``&``) operator.

        Returns
        -------
        response : `~sunpy.net.jsoc.jsoc.JSOCResponse` object
            A collection of records that the query returns.

        Examples
        --------

        *Example 1*

        Request all AIA 304 image data between 2014-01-01T00:00 and
        2014-01-01T01:00::

            >>> import astropy.units as u
            >>> from sunpy.net import jsoc
            >>> from sunpy.net import attrs as a
            >>> client = jsoc.JSOCClient()  # doctest: +REMOTE_DATA
            >>> response = client.search(a.jsoc.Time('2017-09-06T12:00:00', '2017-09-06T12:02:00'),
            ...                          a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Wavelength(304*u.AA),
            ...                          a.jsoc.Segment('image'))  # doctest: +REMOTE_DATA
            >>> print(response)  # doctest: +REMOTE_DATA
                   T_REC         TELESCOP INSTRUME WAVELNTH CAR_ROT
            -------------------- -------- -------- -------- -------
            2017-09-06T11:59:59Z  SDO/AIA    AIA_4      304    2194
            2017-09-06T12:00:11Z  SDO/AIA    AIA_4      304    2194
            2017-09-06T12:00:23Z  SDO/AIA    AIA_4      304    2194
            2017-09-06T12:00:35Z  SDO/AIA    AIA_4      304    2194
            2017-09-06T12:00:47Z  SDO/AIA    AIA_4      304    2194
            2017-09-06T12:00:59Z  SDO/AIA    AIA_4      304    2194
            2017-09-06T12:01:11Z  SDO/AIA    AIA_4      304    2194
            2017-09-06T12:01:23Z  SDO/AIA    AIA_4      304    2194
            2017-09-06T12:01:35Z  SDO/AIA    AIA_4      304    2194
            2017-09-06T12:01:47Z  SDO/AIA    AIA_4      304    2194
            2017-09-06T12:01:59Z  SDO/AIA    AIA_4      304    2194

        *Example 2*

        Request keyword data of ``hmi.v_45s`` for certain specific keywords only::

            >>> import astropy.units as u
            >>> from sunpy.net import jsoc
            >>> from sunpy.net import attrs as a
            >>> client = jsoc.JSOCClient()  # doctest: +REMOTE_DATA
            >>> response = client.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T00:10:00'),
            ...                          a.jsoc.Series('hmi.v_45s'),
            ...                          a.jsoc.Keys('T_REC, DATAMEAN, OBS_VR'))  # doctest: +REMOTE_DATA
            >>> print(response)  # doctest: +REMOTE_DATA
                         T_REC               DATAMEAN            OBS_VR
            ----------------------- ------------------ ------------------
            2014.01.01_00:00:45_TAI        1906.518188        1911.202614
            2014.01.01_00:01:30_TAI        1908.876221        1913.945512
            2014.01.01_00:02:15_TAI          1911.7771 1916.6679989999998
            2014.01.01_00:03:00_TAI        1913.422485 1919.3699239999999
            2014.01.01_00:03:45_TAI        1916.500488        1922.050862
            2014.01.01_00:04:30_TAI        1920.414795 1924.7110050000001
            2014.01.01_00:05:15_TAI        1922.636963         1927.35015
            2014.01.01_00:06:00_TAI 1924.6973879999998        1929.968523
            2014.01.01_00:06:45_TAI        1927.758301 1932.5664510000001
            2014.01.01_00:07:30_TAI        1929.646118         1935.14288
            2014.01.01_00:08:15_TAI        1932.097046        1937.698521
            2014.01.01_00:09:00_TAI 1935.7286379999998         1940.23353
            2014.01.01_00:09:45_TAI        1937.754028        1942.747605
            2014.01.01_00:10:30_TAI 1940.1462399999998        1945.241147

            *Example 3*

        Request data of ``aia.lev1_euv_12s`` on the basis of PrimeKeys other than ``T_REC``::

            >>> import astropy.units as u
            >>> from sunpy.net import jsoc
            >>> from sunpy.net import attrs as a
            >>> client = jsoc.JSOCClient()  # doctest: +REMOTE_DATA
            >>> response = client.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T00:01:00'),
            ...                          a.jsoc.Series('aia.lev1_euv_12s'),
            ...                          a.jsoc.PrimeKey('WAVELNTH','171'))  # doctest: +REMOTE_DATA
            >>> print(response)  # doctest: +REMOTE_DATA
                   T_REC         TELESCOP INSTRUME WAVELNTH CAR_ROT
            -------------------- -------- -------- -------- -------
            2014-01-01T00:00:01Z  SDO/AIA    AIA_3      171    2145
            2014-01-01T00:00:13Z  SDO/AIA    AIA_3      171    2145
            2014-01-01T00:00:25Z  SDO/AIA    AIA_3      171    2145
            2014-01-01T00:00:37Z  SDO/AIA    AIA_3      171    2145
            2014-01-01T00:00:49Z  SDO/AIA    AIA_3      171    2145
            2014-01-01T00:01:01Z  SDO/AIA    AIA_3      171    2145

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
        Builds a jsoc query, similar to query method, and takes similar inputs.

        Complex queries to be easily formed using logical operators such as
        ``&`` and ``|``, in the same way as the query function.

        Parameters
        ----------
        query : a variable number of `~sunpy.net.jsoc.attrs`
                as parameters, which are chained together using
                the ``AND`` (``&``) operator.

        Returns
        -------
        res : `~pandas.DataFrame` object
            A collection of metadata of all the files.

        Example
        -------

        Request metadata or all all AIA 304 image data between 2014-01-01T00:00 and
        2014-01-01T01:00.

        Since, the function only performs a lookdata, and does not make a proper export
        request, attributes like Segment need not be passed::

            >>> import astropy.units as u
            >>> from sunpy.net import jsoc
            >>> from sunpy.net import attrs as a
            >>> client = jsoc.JSOCClient()  # doctest: +REMOTE_DATA
            >>> metadata = client.search_metadata(
            ...                         a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T00:01:00'),
            ...                         a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Wavelength(304*u.AA))  # doctest: +REMOTE_DATA
            >>> print(metadata[['T_OBS', 'WAVELNTH']])  # doctest: +REMOTE_DATA
                                                                        T_OBS  WAVELNTH
            aia.lev1_euv_12s[2014-01-01T00:00:01Z][304]  2014-01-01T00:00:08.57Z       304
            aia.lev1_euv_12s[2014-01-01T00:00:13Z][304]  2014-01-01T00:00:20.58Z       304
            aia.lev1_euv_12s[2014-01-01T00:00:25Z][304]  2014-01-01T00:00:32.57Z       304
            aia.lev1_euv_12s[2014-01-01T00:00:37Z][304]  2014-01-01T00:00:44.58Z       304
            aia.lev1_euv_12s[2014-01-01T00:00:49Z][304]  2014-01-01T00:00:56.57Z       304
            aia.lev1_euv_12s[2014-01-01T00:01:01Z][304]  2014-01-01T00:01:08.59Z       304

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
        Request that JSOC stages the data for download. This method will not
        wait for the request to be staged.

        Parameters
        ----------
        jsoc_response : `~sunpy.net.jsoc.jsoc.JSOCResponse` object
            The results of a query

        Returns
        -------
        requests : `~drms.ExportRequest` object or
                   a list of  `~drms.ExportRequest` objects

            Request Id can be accessed by requests.id
            Request status can be accessed by requests.status

        """

        requests = []
        self.query_args = jsoc_response.query_args
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

            requests.append(r)

        if len(requests) == 1:
            return requests[0]
        return requests

    @deprecated('0.9', alternative='drms.ExportRequest.status')
    def check_request(self, requests):
        """
        Check the status of a request and print out a message about it.

        Parameters
        ----------
        requests :  `~drms.ExportRequest` object or
                   a list of  `~drms.ExportRequest` objects,
                   returned by `~sunpy.net.jsoc.jsoc.JSOCClient.request_data`

        Returns
        -------
        status : `int` or `list`
            A status or list of status' that were returned by JSOC.

        """
        # Convert IDs to a list if not already
        if not isiterable(requests) or isinstance(requests, drms.ExportRequest):
            requests = [requests]

        allstatus = []
        for request in requests:
            status = request.status

            if status == request._status_code_ok:  # Data ready to download
                print("Request {0} was exported at {1} and is ready to "
                      "download.".format(request.id,
                                         request._d['exptime']))
            elif status in request._status_codes_pending:
                print_message = "Request {0} was submitted {1} seconds ago, "\
                                "it is not ready to download."
                print(print_message.format(request.id,
                                           request._d['wait']))
            else:
                print_message = "Request returned status: {0} with error: {1}"
                json_status = request.status
                json_error = request._d.get('error', '')
                print(print_message.format(json_status, json_error))

            allstatus.append(status)

        if len(allstatus) == 1:
            return allstatus[0]
        return allstatus

    def fetch(self, jsoc_response, path=None, overwrite=False, progress=True,
              max_conn=5, downloader=None, sleep=10):
        """
        Make the request for the data in a JSOC response and wait for it to be
        staged and then download the data.

        Parameters
        ----------
        jsoc_response : `~sunpy.net.jsoc.jsoc.JSOCResponse` object
            A response object

        path : `str`
            Path to save data to, defaults to SunPy download dir

        overwrite : `bool`
            Replace files with the same name if True

        progress : `bool`
            Print progress info to terminal

        max_conns : `int`
            Maximum number of download connections.

        downloader: `~sunpy.net.download.Downloader` instance
            A Custom downloader to use

        sleep : `int`
            The number of seconds to wait between calls to JSOC to check the status
            of the request.

        Returns
        -------
        results : a `~sunpy.net.download.Results` instance
            A Results object

        """

        # Make staging request to JSOC
        responses = self.request_data(jsoc_response)
        # Make response iterable
        if not isiterable(responses):
            responses = [responses]
        # Add them to the response for good measure
        jsoc_response.requests = [r for r in responses]
        time.sleep(sleep/2.)

        r = Results(lambda x: None, done=lambda maps: [v['path'] for v in maps.values()])

        for response in responses:
            response.wait(verbose=progress)
            r = self.get_request(response, path=path, overwrite=overwrite,
                                 progress=progress, results=r)

        return r

    @deprecated('0.8', alternative='JSOCClient.fetch')
    def get(self, jsoc_response, path=None, overwrite=False, progress=True,
            max_conn=5, downloader=None, sleep=10):
        """
        See `~sunpy.net.jsoc.jsoc.JSOCClient.fetch`
        """
        return self.fetch(jsoc_response, path=path, overwrite=overwrite, progress=progress,
                          max_conn=max_conn, downloader=downloader, sleep=sleep)

    def get_request(self, requests, path=None, overwrite=False, progress=True,
                    max_conn=5, downloader=None, results=None):
        """
        Query JSOC to see if the request(s) is ready for download.

        If the request is ready for download, it will then download it.

        Parameters
        ----------
        requests : `~drms.ExportRequest`, `str`, `list`
            `~drms.ExportRequest` objects or `str` request IDs or lists
            returned by `~sunpy.net.jsoc.jsoc.JSOCClient.request_data`.

        path : `str`
            Path to save data to, defaults to SunPy download dir.

        overwrite : `bool`
            Replace files with the same name if True.

        progress : `bool`
            Print progress info to terminal.

        max_conns : `int`
            Maximum number of download connections.

        downloader : `~sunpy.net.download.Downloader`
            A Custom downloader to use

        results: `~sunpy.net.download.Results`
            A `~sunpy.net.download.Results` manager to use.

        Returns
        -------
        res: `~sunpy.net.download.Results`
            A `~sunpy.net.download.Results` instance or `None` if no URLs to download

        """
        c = drms.Client()

        # Convert Responses to a list if not already
        if isinstance(requests, six.string_types) or not isiterable(requests):
            requests = [requests]

        # Ensure all the requests are drms ExportRequest objects
        for i, request in enumerate(requests):
            if isinstance(request, six.string_types):
                r = c.export_from_id(request)
                requests[i] = r

        # We only download if all are finished
        if not all([r.has_succeeded() for r in requests]):
            raise NotExportedError("Can not download as not all the requests"
                                   "have been exported for download yet.")

        # Ensure path has a {file} in it
        if path is None:
            default_dir = config.get("downloads", "download_dir")
            path = os.path.join(default_dir, '{file}')
        elif isinstance(path, six.string_types) and '{file}' not in path:
            path = os.path.join(path, '{file}')

        paths = []
        for request in requests:
            for filename in request.data['filename']:
                # Ensure we don't duplicate the file extension
                ext = os.path.splitext(filename)[1]
                if path.endswith(ext):
                    fname = path.strip(ext)
                else:
                    fname = path
                fname = fname.format(file=filename)
                fname = os.path.expanduser(fname)
                fname = partial(simple_path, fname)
                paths.append(fname)

        if downloader is None:
            downloader = Downloader(max_conn=max_conn, max_total=max_conn)

        # A Results object tracks the number of downloads requested and the
        # number that have been completed.
        if results is None:
            results = Results(lambda _: downloader.stop(),
                              done=lambda maps: [v['path'] for v in maps.values()])

        urls = []
        for request in requests:

            if request.status == 0:
                for index, data in request.data.iterrows():
                    is_file = os.path.isfile(paths[index].args[0])
                    if overwrite or not is_file:
                        url_dir = request.request_url + '/'
                        urls.append(urllib.parse.urljoin(url_dir, data['filename']))

                    if not overwrite and is_file:
                        print_message = "Skipping download of file {} as it " \
                                        "has already been downloaded. " \
                                        "If you want to redownload the data, "\
                                        "please set overwrite to True"
                        print(print_message.format(data['filename']))
                        # Add the file on disk to the output
                        results.map_.update({data['filename']:
                                            {'path': paths[index].args[0]}})
        if urls:
            if progress:
                print_message = "{0} URLs found for download. Full request totalling {1}MB"
                print(print_message.format(len(urls), request._d['size']))
            for i, url in enumerate(urls):
                downloader.download(url, callback=results.require([url]),
                                    errback=lambda x: print(x), path=paths[i])

        else:
            # Make Results think it has finished.
            results.require([])
            results.poke()

        return results

    def _make_recordset(self, series, start_time='', end_time='', wavelength='',
                        segment='', primekey={}, **kwargs):
        """
        Take the query arguments and build a record string.

        All the primekeys are now stored in primekey dict, including Time and Wavelength
        which were passed through pre-defined attributes. The following piece of code,
        extracts the passed prime-keys and arranges it in the order as it appears in the
        JSOC database.

        `pkeys_isTime` is a Pandas DataFrame, whose index values are the Prime-key names
        and the column stores a boolean value, identifying whether the prime-key is a
        Time-type prime-key or not. Since, time-type prime-keys exist by different names,
        we made it uniform in the above piece of code, by storing the time-type primekey
        with a single name `TIME`.

        Considering an example, if the primekeys that exist for a given series are
        ['HARPNUM', 'T_OBS', 'WAVELNTH'], we will consider three different cases of the
        passed primekeys.

        pkeys_isTime.index.values = ['HARPNUM', 'T_OBS', 'WAVELNTH']

        Case 1
        ------

        primekey = {'T_OBS' : , '2014.01.01_00:00:45_TAI',
                    'HARPNUM' : '4864',
                    'WAVELNTH': '605'}

        If the primekey dict is as above, then pkstr should be as:

        pkstr = '{4864}{2014.01.01_00:00:45_TAI}{605}'

        Case 2
        ------

        primekey = {'T_OBS' : , '2014.01.01_00:00:45_TAI',
                    'WAVELNTH': '605'}

        If the primekey dict is as above, then pkstr should be as:

        pkstr = '{}{2014.01.01_00:00:45_TAI}{605}'

        Case 3
        ------

        primekey = {'T_OBS' : , '2014.01.01_00:00:45_TAI'}

        If the primekey dict is as above, then pkstr should be as:

        pkstr = '{}{2014.01.01_00:00:45_TAI}'

        The idea behind this should be clear. We build up the `pkstr` string
        containing the values of the prime-keys passed in the same order as
        it occurs in the list `pkeys_isTime.index.values`, i.e. how it is stored
        in the online database. Any missing prime-keys should be compensated by
        an empty {}, if it occurs before any passed prime-key. Any empty curly braces
        that is present at last of the pkstr, can be skipped.

        """

        # Extract and format segment
        # Convert list of segments into a comma-separated string
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
            # Check whether any primekey listed in PKEY_LIST_TIME has been passed through
            # PrimeKey() attribute. If yes, raise an error, since Time can only be passed
            # either through PrimeKey() attribute or Time() attribute.
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
            # This is executed when Time has not been passed through Time() attribute.
            # `match` stores all the time-type prime-keys that has been passed through
            # PrimeKey() attribute. The length of `match` won't ever be greater than 1,
            # but it is a good idea to keep a check.
            match = set(primekey.keys()) & PKEY_LIST_TIME
            if len(match) > 1:
                error_message = "Querying of series, having more than 1 Time-type "\
                                "prime-keys is not yet supported. Alternative is to "\
                                "use only one of the primekey to query for series data."
                raise ValueError(error_message)

            if match:
                timestr = '{0}'.format(primekey.pop(list(match)[0], ''))
            else:
                timestr = ''

        if wavelength:
            if not primekey.get('WAVELNTH', ''):
                if isinstance(wavelength, list):
                    wavelength = [int(np.ceil(wave.to(u.AA).value)) for wave in wavelength]
                    wavelength = str(wavelength)
                else:
                    wavelength = '{0}'.format(int(np.ceil(wavelength.to(u.AA).value)))

            else:
                # This is executed when wavelength has been passed both through PrimeKey()
                # and Wavelength().
                error_message = "Wavelength attribute has been passed both as a Wavelength()"\
                                " and PrimeKey(). Please provide any one of them"\
                                " or separate them by OR operator."
                raise ValueError(error_message)

        else:
            # This is executed when wavelength has been passed through PrimeKey().
            wavelength = '{0}'.format(primekey.pop('WAVELNTH', ''))

        # Populate primekey dict with formatted Time and Wavlength.
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
            # The loop is iterating over the list of prime-keys existing for the given series.
            if len(primekey) > 0:
                if pkeys_isTime[pkey]:
                    pkstr += '[{0}]'.format(primekey.pop('TIME', ''))
                else:
                    pkstr += '[{0}]'.format(primekey.pop(pkey, ''))
            else:
                break
                # break because we can skip adding {} at the end of pkstr, if the primekey
                # dict is empty.

        if not pkstr:
            # pkstr cannot be totally empty
            error_message = "Atleast one PrimeKey must be passed."
            raise ValueError(error_message)

        dataset = '{series}{primekey}{segment}'.format(series=series,
                                                       primekey=pkstr,
                                                       segment=segment)

        return dataset

    def _lookup_records(self, iargs):
        """
        Do a LookData request to JSOC to workout what results the query returns.
        """

        keywords_default = ['T_REC', 'TELESCOP', 'INSTRUME', 'WAVELNTH', 'CAR_ROT']
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
        # Get a set of the PrimeKeys that exist for the given series, and check
        # whether the passed PrimeKeys is a subset of that.
        pkeys = c.pkeys(iargs['series'])
        pkeys_passed = iargs.get('primekey', None)  # pkeys_passes is a dict, with key-value pairs.
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
        # Get a set of the segments that exist for the given series, and check
        # whether the passed segments is a subset of that.
        si = c.info(iargs['series'])
        segs = list(si.segments.index.values)          # Fetches all valid segment names
        segs_passed = iargs.get('segment', None)
        if segs_passed is not None:

            if not isinstance(segs_passed, list) and not isinstance(segs_passed, six.string_types):
                error_message = "Segments can only be passed as a comma-separated"\
                                " string or a list of strings."
                raise TypeError(error_message)

            elif isinstance(segs_passed, six.string_types):
                segs_passed = segs_passed.replace(' ', '').split(',')

            if not set(segs_passed) <= set(segs):
                error_message = "Unexpected Segments were passed. The series {series} "\
                                "contains the following Segments {segs}"
                raise ValueError(error_message.format(series=iargs['series'], segs=segs))

            iargs['segment'] = segs_passed

        # If Time has been passed as a PrimeKey, convert the Time object into TAI time scale,
        # and then, convert it to datetime object.

        iargs['start_time'] = iargs['start_time'].tai.datetime
        iargs['end_time'] = iargs['end_time'].tai.datetime

        ds = self._make_recordset(**iargs)

        # Convert the list of keywords into comma-separated string.
        if isinstance(keywords, list):
            key = str(keywords)[1:-1].replace(' ', '').replace("'", '')
        else:
            key = keywords

        r = c.query(ds, key=key, rec_index=isMeta)

        # If the method was called from search_metadata(), return a Pandas Dataframe,
        # otherwise return astropy.table
        if isMeta:
            return r

        if r is None or r.empty:
            return astropy.table.Table()
        else:
            return astropy.table.Table.from_pandas(r)

    @classmethod
    def _can_handle_query(cls, *query):
        chkattr = ['Series', 'Protocol', 'Notify', 'Wavelength', 'Time',
                   'Segment', 'Keys', 'PrimeKey', 'Sample']

        return all([x.__class__.__name__ in chkattr for x in query])
