# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 13:44:56 2014

@author: stuart
"""
from __future__ import print_function, absolute_import

import os
import time
import urlparse
import warnings

import requests
import astropy.time

from sunpy import config
from sunpy.time import parse_time, TimeRange
from sunpy.net.download import Downloader
from sunpy.net.vso.vso import Results

__all__ = ['JSOCClient']

JSOC_URL = 'http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_fetch'
BASE_DL_URL = 'http://jsoc.stanford.edu'

class JSOCClient(object):
    """
    This is a Client to the JSOC Data Export service.

    This is not warm and fluffy like the VSO, this is hardcore crazy. It is
    more efficient for large scale SDO queries than the VSO, however harder to
    use.

    Notes
    -----
    This Client mocks input to this site: http://jsoc.stanford.edu/ajax/exportdata.html
    Therefore that is a good resource if things are mis-behaving.
    The full list of 'series' is availible through this site: http://jsoc.stanford.edu/

    You can build more complex queries by specifiying parameters to POST to JSOC via keyword
    arguments. You can generate these kwargs using the Export Data page at JSOC.

    Examples
    --------
    Query JSOC for some HMI data at 45 second cadence

    >>> from sunpy.net import jsoc
    >>> client = jsoc.JSOCClient()
    >>> IDs = client.jsoc_query('2012/1/1T00:00:00', '2012/1/1T00:00:45', 'hmi.m_45s')

    The returned `IDs` is a list of JSOC request identifiers, you can check the status of
    a request thus:

    >>> status = client.check_status(IDs)

    Once the request has been staged (Status 1) you can download the data:

    >>> res = client.get(IDs)

    This returns a Results instance which can be used to watch the progress
    of the download.
    """

    def jsoc_query(self, start_time, end_time, series, **kwargs):
        """
        Build a JSOC query and submit it to JSOC for processing.

        Parameters
        ----------
        start_time: datetime, astropy.time or string
            The time in UTC (or astropy.time object) specifying the start time.

        end_time: datetime, astropy.time or string
            The time in UTC (or astropy.time object) specifying the end time.

        series: string
            A string representing the JSOC data series to download e.g. 'hmi.M_45s'

        notify: string
            An email address to get a notification to when JSOC has staged your request

        protocol: string
            The type of download to request one of ("FITS", "JPEG", "MPG", "MP4", or "as-is").
            Only FITS is supported, the others will require extra keywords.

        compression: string
            'rice' or None, download FITS files with RICE compression.

        kwargs: dict
            Extra keywords to put into the POST payload to JSOC.

        Returns
        -------
        requestIDs: list
            A list of the requestIDs generated from your query
        """

        # A little (hidden) debug feature
        return_resp = kwargs.pop('return_resp', False)

        start_time = self._process_time(start_time)
        end_time = self._process_time(end_time)

        # Do a multi-request
        responses = self._multi_request(start_time, end_time, series, **kwargs)

        for i, response in enumerate(responses):
            #TODD: catch non 200 return
            if response.json()['status'] != 2:
                warnings.warn(
                Warning("Query {0} retuned status {1} with error {2}".format(i,
                                                     response.json()['status'],
                                                     response.json()['error'])))
                responses.pop(i)
        #Extract the IDs from the JSON
        requestIDs = [response.json()['requestid'] for response in responses]

        if return_resp:
            return responses

        return requestIDs

    def check_request(self, requestIDs):
        """
        Check the status of a request and print out a messgae about it

        Parameters
        ----------
        requestIDs: list or string
            A list of requestIDs to check

        Returns
        -------
        status: list
            A list of status' that were returned by JSOC
        """
        # Convert IDs to a list if not already
        if not astropy.utils.misc.isiterable(requestIDs) or isinstance(requestIDs, basestring):
            requestIDs = [requestIDs]

        allstatus = []
        for request_id in requestIDs:
            u = self._request_status(request_id)
            status = int(u.json()['status'])

            if status == 0: #Data ready to download
                print("Request {0} was exported at {1} and is ready to download.".format(u.json()['requestid'],
                                                                                       u.json()['exptime']))
            elif status == 1:
                print("Request {0} was submitted {1} seconds ago, it is not ready to download.".format(
                                                             u.json()['requestid'], u.json()['wait']))
            else:
                print("Request returned status: {0} with error: {1}".format(
                                    u.json()['status'], u.json()['error']))

            allstatus.append(status)

        return allstatus


    def wait_get(self, requestIDs, path=None, overwrite=False, progress=True,
            max_conn=5, sleep=10):
        """
        Same as get() excepts it will wait until the download has been staged.

        Parameters
        ----------
        requestIDs: list or string
            One or many requestID strings

        path: string
            Path to save data to, defaults to SunPy download dir

        overwrite: bool
            Replace files with the same name if True

        progress: bool
            Print progress info to terminal

        max_conns: int
            Maximum number of download connections.

        downloader: sunpy.download.Downloder instance
            A Custom downloader to use

        sleep: int
            The number of seconds to wait between calls to JSOC to check the status
            of the request.

        Returns
        -------
        downloader: a sunpy.net.download.Downloader instance
            A Downloader instance
        """
        # Convert IDs to a list if not already
        if not astropy.utils.misc.isiterable(requestIDs) or isinstance(requestIDs, basestring):
            requestIDs = [requestIDs]

        r = Results(lambda x: None)

        while requestIDs:
            for i, request_id in enumerate(requestIDs):
                u = self._request_status(request_id)

                if progress:
                    self.check_request(request_id)

                if u.status_code == 200 and u.json()['status'] == '0':
                    rID = requestIDs.pop(i)
                    r = self.get(rID, path=path, overwrite=overwrite,
                             progress=progress, results=r)

                else:
                    time.sleep(sleep)

        return r

    def get(self, requestIDs, path=None, overwrite=False, progress=True,
            max_conn=5, downloader=None, results=None):
        """
        Query JSOC to see if request_id is ready for download.

        If the request is ready for download download it.

        Parameters
        ----------
        requestIDs: list or string
            One or many requestID strings

        path: string
            Path to save data to, defaults to SunPy download dir

        overwrite: bool
            Replace files with the same name if True

        progress: bool
            Print progress info to terminal

        max_conns: int
            Maximum number of download connections.

        downloader: sunpy.download.Downloder instance
            A Custom downloader to use

        results: Results instance
            A Results manager to use.

        Returns
        -------
        res: Results
            A Results instance or None if no URLs to download
        """

        # Convert IDs to a list if not already
        if not astropy.utils.misc.isiterable(requestIDs) or isinstance(requestIDs, basestring):
            requestIDs = [requestIDs]

        if path is None:
            path = config.get('downloads','download_dir')

        if downloader is None:
            downloader = Downloader(max_conn=max_conn, max_total=max_conn)

        # A Results object tracks the number of downloads requested and the
        # number that have been completed.
        if results is None:
            results = Results(lambda x: None)

        urls = []
        for request_id in requestIDs:
            u = self._request_status(request_id)

            if u.status_code == 200 and u.json()['status'] == '0':
                for ar in u.json()['data']:
                    if overwrite or not os.path.isfile(os.path.join(path, ar['filename'])):
                        urls.append(urlparse.urljoin(BASE_DL_URL + u.json()['dir']+'/', ar['filename']))
                if progress:
                    print("{0} URLs found for Download. Totalling {1}MB".format(len(urls), u.json()['size']))

            else:
                if progress:
                    self.check_request(request_id)

        if urls:
            for url, rcall in list(zip(urls, list(map(lambda x: results.require([x]), urls)))):
                downloader.download(url, callback=rcall, path=path)

        else:
            #Make Results think it has finished.
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
        time = astropy.time.Time(time, scale='tai')

        return time.datetime


    def _make_query_payload(self, start_time, end_time, series, notify='',
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

        payload = {'ds':'{0}[{1}-{2}]'.format(series, start_time.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                                           end_time.strftime("%Y.%m.%d_%H:%M:%S_TAI")),
                   'format':'json',
                   'method':'url',
                   'notify':notify,
                   'op':'exp_request',
                   'process':'n=0|no_op',
                   'protocol':jprotocol,
                   'requestor':'none',
                   'filenamefmt':'{0}.{{T_REC:A}}.{{CAMERA}}.{{segment}}'.format(series)}

        payload.update(kwargs)

        return payload

    def _send_jsoc_request(self, start_time, end_time, series, notify='',
                          protocol='FITS', compression='rice', **kwargs):
        """
        Request that JSOC stages data for download

        This routine puts in a POST request to JSOC
        """

        payload = self._make_query_payload(start_time, end_time, series, notify=notify,
                          protocol=protocol, compression=compression, **kwargs)

        r = requests.post(JSOC_URL, data=payload)

        if r.status_code != 200:
            raise Exception("JSOC POST Request returned code {0}".format(r.status_code))

        return r, r.json()

    def _multi_request(self, start_time, end_time, series, **kwargs):
        """
        Make a series of requests to avoid the 100GB limit
        """
        tr = TimeRange(start_time, end_time)
        returns = []

        response, json_response = self._send_jsoc_request(start_time, end_time, series, **kwargs)

        if json_response['status'] == 3 and json_response['error'] == 'Request exceeds max byte limit of 100000MB':
            returns.append(self._multi_request(tr.start(), tr.center(), series, **kwargs)[0])
            returns.append(self._multi_request(tr.center(), tr.end(), series, **kwargs)[0])
        else:
            returns.append(response)

        return returns

    def _request_status(self, request_id):
        """
        GET the status of a request ID
        """
        payload = {'op':'exp_status', 'requestid':request_id}
        u = requests.get(JSOC_URL, params=payload)

        return u
