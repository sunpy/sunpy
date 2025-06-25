import os
import copy
import json
import time
import urllib
from pathlib import Path

import drms
import numpy as np
import pandas as pd
import parfive
from packaging.version import Version

import astropy.table
import astropy.time
import astropy.units as u

from sunpy import config, log
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable, convert_row_to_table
from sunpy.net.jsoc.attrs import walker
from sunpy.util.exceptions import warn_user
from sunpy.util.parfive_helpers import Downloader, Results

__all__ = ['JSOCClient', 'JSOCResponse']


PRIMEKEY_LIST_TIME = {'T_START', 'T_REC', 'T_OBS', 'MidTime', 'OBS_DATE',
                      'obsdate', 'DATE_OBS', 'starttime', 'stoptime', 'UTC_StartTime'}
parfive_version = Version(parfive.__version__)


class NotExportedError(Exception):
    pass


class JSOCResponse(QueryResponseTable):
    query_args = astropy.table.TableAttribute()
    requests = astropy.table.TableAttribute()
    display_keys = ['T_REC', 'TELESCOP', 'INSTRUME', 'WAVELNTH', 'CAR_ROT']
    # This variable is used to detect if the result has been sliced before it is passed
    # to fetch and issue a warning to the user about not being able to post-filter JSOC searches.
    _original_num_rows = astropy.table.TableAttribute(default=None)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._original_num_rows = len(self)


class JSOCClient(BaseClient):
    """
    Provides access to the JSOC Data Export service.

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
    The full list of ``Series`` is available via ` this page <http://jsoc.stanford.edu/JsocSeries_DataProducts_map.html>`__.

    JSOC requires a validated email address, you can pass in your validated email address
    using the `~sunpy.net.jsoc.attrs.Notify` attribute. You have to register your email address
    with JSOC beforehand `here <http://jsoc.stanford.edu/ajax/register_email.html>`__.

    The backend of SunPy's JSOC Client uses `drms package <https://github.com/sunpy/drms>`__.
    The tutorials can be :ref:`found here <sunpy-tutorial-acquiring-data-jsoc>`.
    This can be used to build complex queries, by directly inputting the query string.

    Examples
    --------

    *Example 1*

    Query JSOC for some HMI data at 45 second cadence::

        >>> from sunpy.net import jsoc
        >>> from sunpy.net import attrs as a
        >>> client = jsoc.JSOCClient()
        >>> response = client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T00:10:00'),
        ...                          a.jsoc.Series('hmi.m_45s'), a.jsoc.Notify("sunpy@sunpy.org"))  # doctest: +SKIP

        The response object holds the records that your query will return:

        >>> print(response)   # doctest: +SKIP
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

        >>> client = jsoc.JSOCClient()  # doctest: +SKIP
        >>> response = client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T00:10:00'),
        ...                          a.jsoc.Series('hmi.m_45s'))  # doctest: +SKIP

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
        >>> client = jsoc.JSOCClient()  # doctest: +SKIP
        >>> response = client.search(a.Time('2014/1/1T00:00:00', '2014/1/1T00:00:36'),
        ...                          a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Segment('image'),
        ...                          a.Wavelength(171*u.AA), a.jsoc.Notify("sunpy@sunpy.org"))  # doctest: +SKIP

        The response object holds the records that your query will return:

        >>> print(response)  # doctest: +SKIP
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
        2


    The status codes represent the following:

    0: Completed, ready to download.
    2: Pending, but downloadable.
    6: Not found.

    Once the status code is 0 or 2, you can download the data using the ``get_request`` method::

        >>> res = client.get_request(requests)  # doctest: +SKIP

    This returns a Results instance which can be used to watch the progress
    of the download::

        >>> res.wait(progress=True)   # doctest: +SKIP

    """
    # Default number of max connections that the Downloader opens
    default_max_conn = 1

    @property
    def info_url(self):
        return 'http://jsoc.stanford.edu'

    def search(self, *query, **kwargs):
        """
        Build a JSOC query and submit it to JSOC for processing.

        Takes a variable number of `~sunpy.net.jsoc.attrs` as parameters,
        which are chained together using the AND (``&``) operator.

        Complex queries to be easily formed using logical operators such as
        ``&`` and ``|``, in the same way as the VSO client.

        Parameters
        ----------
        *query : a variable number of `~sunpy.net.jsoc.attrs`
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
            >>> client = jsoc.JSOCClient()  # doctest: +SKIP
            >>> response = client.search(a.Time('2017-09-06T12:00:00', '2017-09-06T12:02:00'),
            ...                          a.jsoc.Series('aia.lev1_euv_12s'), a.Wavelength(304*u.AA),
            ...                          a.jsoc.Segment('image'))  # doctest: +SKIP
            >>> print(response)  # doctest: +SKIP
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

        Request keyword data of "hmi.v_45s" and show specific columns only::

            >>> import astropy.units as u
            >>> from sunpy.net import jsoc
            >>> from sunpy.net import attrs as a
            >>> client = jsoc.JSOCClient()  # doctest: +SKIP
            >>> response = client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T00:10:00'),
            ...                          a.jsoc.Series('hmi.v_45s'))  # doctest: +SKIP
            >>> print(response.show('T_REC', 'WAVELNTH', 'CAR_ROT'))  # doctest: +SKIP
                     T_REC          WAVELNTH CAR_ROT
            ----------------------- -------- -------
            2014.01.01_00:00:45_TAI   6173.0    2145
            2014.01.01_00:01:30_TAI   6173.0    2145
            2014.01.01_00:02:15_TAI   6173.0    2145
            2014.01.01_00:03:00_TAI   6173.0    2145
            2014.01.01_00:03:45_TAI   6173.0    2145
            2014.01.01_00:04:30_TAI   6173.0    2145
            2014.01.01_00:05:15_TAI   6173.0    2145
            2014.01.01_00:06:00_TAI   6173.0    2145
            2014.01.01_00:06:45_TAI   6173.0    2145
            2014.01.01_00:07:30_TAI   6173.0    2145
            2014.01.01_00:08:15_TAI   6173.0    2145
            2014.01.01_00:09:00_TAI   6173.0    2145
            2014.01.01_00:09:45_TAI   6173.0    2145
            2014.01.01_00:10:30_TAI   6173.0    2145

        *Example 3*

        Request data of ``aia.lev1_euv_12s`` on the basis of PrimeKeys other than ``T_REC``::

            >>> import astropy.units as u
            >>> from sunpy.net import jsoc
            >>> from sunpy.net import attrs as a
            >>> client = jsoc.JSOCClient()  # doctest: +SKIP
            >>> response = client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T00:01:00'),
            ...                          a.jsoc.Series('aia.lev1_euv_12s'),
            ...                          a.jsoc.PrimeKey('WAVELNTH','171'))  # doctest: +SKIP
            >>> print(response)  # doctest: +SKIP
                   T_REC         TELESCOP INSTRUME WAVELNTH CAR_ROT
            -------------------- -------- -------- -------- -------
            2014-01-01T00:00:01Z  SDO/AIA    AIA_3      171    2145
            2014-01-01T00:00:13Z  SDO/AIA    AIA_3      171    2145
            2014-01-01T00:00:25Z  SDO/AIA    AIA_3      171    2145
            2014-01-01T00:00:37Z  SDO/AIA    AIA_3      171    2145
            2014-01-01T00:00:49Z  SDO/AIA    AIA_3      171    2145
            2014-01-01T00:01:01Z  SDO/AIA    AIA_3      171    2145

        """
        return_results = JSOCResponse(client=self)
        query = and_(*query)
        blocks = []
        for block in walker.create(query):
            extra_kwargs = kwargs.copy()
            extra_kwargs |= block
            # Update blocks with deep copy of extra_kwargs because in _make_recordset we use .pop() on
            # element from extra_kwargs
            blocks.append(copy.deepcopy(extra_kwargs))
            return_results = astropy.table.vstack([return_results, self._lookup_records(extra_kwargs)])
        return_results.query_args = blocks
        return_results._original_num_rows = len(return_results)
        return return_results

    def request_data(self, jsoc_response, method='url', **kwargs):
        """
        Request that JSOC stages the data for download. This method will not
        wait for the request to be staged.

        Parameters
        ----------
        jsoc_response : `~sunpy.net.jsoc.jsoc.JSOCResponse` object
            The results of a query
        method : {``"url"``, ``"url-tar"``, ``"url-quick"``}
            Method for requesting JSOC data, can be 'url-tar', 'url' (the default) and 'url-quick'
            If 'url-tar' it will request JSOC to provide single .tar file which contains all data
            If 'url' it will request JSOC to provide all data as separate .fits files
            If 'url-quick' (only with protocol 'as-is') provide all data as separate files,
            but only if data is online.

        Returns
        -------
        requests : `~drms.ExportRequest` object or
                   a list of  `~drms.ExportRequest` objects
            Request Id can be accessed by requests.id
            Request status can be accessed by requests.status
        """

        requests = []
        self.query_args = jsoc_response.query_args
        supported_protocols = sorted(['fits', 'as-is'])
        supported_methods = sorted(['url-tar', 'url', 'url-quick'])
        for block in jsoc_response.query_args:

            ds = self._make_recordset(**block)
            cd = drms.Client(email=block.get('notify', ''))
            protocol = block.get('protocol', 'fits')
            cutout = block.get('cutout')

            if protocol not in supported_protocols:
                error_message = f"Protocols other than {','.join(supported_protocols)} are not supported."
                raise TypeError(error_message)
            if method not in supported_methods:
                error_message = f"Methods other than {','.join(supported_methods)}  are not supported."
                raise TypeError(error_message)
            process = {'im_patch': cutout} if cutout is not None else None

            if method != 'url-tar':
                method = 'url' if protocol == 'fits' else 'url_quick'
            r = cd.export(ds, method=method, protocol=protocol, process=process)

            requests.append(r)

        return requests[0] if len(requests) == 1 else requests

    @convert_row_to_table
    def fetch(self, jsoc_response, path=None, progress=True, overwrite=False,
              downloader=None, wait=True, sleep=10, max_conn=default_max_conn, **kwargs):
        """
        Make the request for the data in a JSOC response and wait for it to be
        staged and then download the data.

        .. note::

            **Only complete searches can be downloaded from JSOC**, this means
            that no slicing operations performed on the results object will
            affect the number of files downloaded.

        Parameters
        ----------
        jsoc_response : `~sunpy.net.jsoc.jsoc.JSOCResponse` object
            A response object
        path : `str`
            Path to save data to, defaults to SunPy download dir
        progress : `bool`, optional
            If `True` show a progress bar showing how many of the total files
            have been downloaded. If `False`, no progress bar will be shown.
        overwrite : `bool` or `str`, optional
            Determine how to handle downloading if a file already exists with the
            same name. If `False` the file download will be skipped and the path
            returned to the existing file, if `True` the file will be downloaded
            and the existing file will be overwritten, if ``'unique'`` the filename
            will be modified to be unique.
        max_conn : `int`
            Maximum number of download connections.
        downloader : `parfive.Downloader`, optional
            The download manager to use.
        wait : `bool`, optional
            If `False` ``downloader.download()`` will not be called. Only has
            any effect if ``downloader`` is not `None`.
        sleep : `int`
            The number of seconds to wait between calls to JSOC to check the status
            of the request.

        Returns
        -------
        results : a `parfive.Results` instance
            A `parfive.Results` object.

        """
        if len(jsoc_response) == 0:
            return parfive.Results()

        for resp in jsoc_response.query_args:
            if 'notify' not in resp:
                raise ValueError('A registered email is required to get data from JSOC. '
                                 'Please supply an email with attrs.jsoc.Notify to Fido.search. '
                                 'Then pass those new results back into Fido.fetch')

        if len(jsoc_response) != jsoc_response._original_num_rows:
            warn_user("Downloading of sliced JSOC results is not supported. "
                      "All the files present in the original response will "
                      "be downloaded when passed to fetch().")

        # Make staging request to JSOC
        responses = self.request_data(jsoc_response)

        defaults = {'max_splits': 1} | kwargs
        # Make response iterable
        if not np.iterable(responses):
            responses = [responses]

        # Add them to the response for good measure
        jsoc_response.requests = list(responses)
        time.sleep(sleep/2.)

        for response in responses:
            response.wait()

        return self.get_request(responses, path=path, overwrite=overwrite,
                                progress=progress, downloader=downloader,
                                wait=wait, max_conn=max_conn, **defaults)

    def get_request(self, requests, path=None, overwrite=False, progress=True,
                    downloader=None, wait=True, max_conn=default_max_conn, **kwargs):
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
        progress : `bool`, optional
            If `True` show a progress bar showing how many of the total files
            have been downloaded. If `False`, no progress bar will be shown.
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
        res: `parfive.Results`
            A `parfive.Results` instance or `None` if no URLs to download
        """
        c = drms.Client()

        # Private communication from JSOC say we should not use more than one connection.
        max_splits = kwargs.get('max_splits', 1)
        if max_splits != 1:
            log.info("Setting max_splits to it's maximum allowed value of 1 for requests made by the JSOCClient.")
            max_splits = 1

        # Convert Responses to a list if not already
        if isinstance(requests, str) or not np.iterable(requests):
            requests = [requests]

        # Ensure all the requests are drms ExportRequest objects
        for i, request in enumerate(requests):
            if isinstance(request, str):
                r = c.export_from_id(request)
                requests[i] = r

        # We only download if all are finished
        if not all([r.has_succeeded() for r in requests]):
            raise NotExportedError("Can not download as not all the requests "
                                   "have been exported for download yet.")

        # Ensure path has a {file} in it
        if path is None:
            default_dir = config.get("downloads", "download_dir")
            path = os.path.join(default_dir, '{file}')
        elif isinstance(path, Path):
            path = str(path)

        if isinstance(path, str) and '{file}' not in path:
            path = os.path.join(path, '{file}')

        paths = []
        for request in requests:
            if request.method == 'url-tar':
                fname = path.format(file=Path(request.tarfile).name)
                paths.append(os.path.expanduser(fname))
            else:
                for filename in request.data['filename']:
                    # Ensure we don't duplicate the file extension
                    ext = os.path.splitext(filename)[1]
                    if path.endswith(ext):
                        fname = path.strip(ext)
                    else:
                        fname = path
                    fname = fname.format(file=filename)
                    fname = os.path.expanduser(fname)
                    paths.append(fname)

        dl_set = True
        if not downloader:
            dl_set = False
            # Private communication from JSOC say we should not use more than one connection.
            if max_conn != self.default_max_conn:
                log.info("Setting max parallel downloads to 1 for the JSOC client.")
            downloader = Downloader(max_conn=max_conn, progress=progress, overwrite=overwrite, max_splits=max_splits)
        urls = []
        for request in requests:
            if request.status == 0:
                if request.protocol == 'as-is' or request.method == 'url-tar':
                    urls.extend(list(request.urls.url))
                else:
                    for _, data in request.data.iterrows():
                        url_dir = request.request_url + '/'
                        urls.append(urllib.parse.urljoin(url_dir, data['filename']))

        if urls:
            if progress:
                log.info(f"{len(urls)} URLs found for download. Full request totaling {request._d['size']}MB")
            for aurl, fname in zip(urls, paths):
                downloader.enqueue_file(aurl, filename=fname, max_splits=max_splits)

        if dl_set and not wait:
            return Results()

        results = downloader.download()
        return results

    def _make_recordset(self, series, start_time='', end_time='', wavelength='',
                        segment='', primekey={}, keyword={}, **kwargs):
        """
        Take the query arguments and build a record string.

        All the primekeys are now stored in primekey dict, including Time and Wavelength
        which were passed through pre-defined attributes. The following piece of code,
        extracts the passed prime-keys and arranges it in the order as it appears in the
        JSOC database.

        `primekeys_istime` is a Pandas DataFrame, whose index values are the Prime-key names
        and the column stores a boolean value, identifying whether the prime-key is a
        Time-type prime-key or not. Since, time-type prime-keys exist by different names,
        we made it uniform in the above piece of code, by storing the time-type primekey
        with a single name `TIME`.

        Considering an example, if the primekeys that exist for a given series are
        ['HARPNUM', 'T_OBS', 'WAVELNTH'], we will consider three different cases of the
        passed primekeys.

        primekeys_istime.index.values = ['HARPNUM', 'T_OBS', 'WAVELNTH']

        Case 1:

        primekey = {'T_OBS' : , '2014.01.01_00:00:45_TAI',
                    'HARPNUM' : '4864',
                    'WAVELNTH': '605'}

        If the primekey dict is as above, then pkstr should be as:

        pkstr = '{4864}{2014.01.01_00:00:45_TAI}{605}'

        Case 2:

        primekey = {'T_OBS' : , '2014.01.01_00:00:45_TAI',
                    'WAVELNTH': '605'}

        If the primekey dict is as above, then pkstr should be as:

        pkstr = '{}{2014.01.01_00:00:45_TAI}{605}'

        Case 3:

        primekey = {'T_OBS' : , '2014.01.01_00:00:45_TAI'}

        If the primekey dict is as above, then pkstr should be as:

        pkstr = '{}{2014.01.01_00:00:45_TAI}'

        The idea behind this should be clear. We build up the `pkstr` string
        containing the values of the prime-keys passed in the same order as
        it occurs in the list `primekeys_istime.index.values`, i.e. how it is stored
        in the online database. Any missing prime-keys should be compensated by
        an empty {}, if it occurs before any passed prime-key. Any empty curly braces
        that is present at last of the pkstr, can be skipped.
        """
        # Extract and format segment
        # Convert list of segments into a comma-separated string
        if segment:
            if isinstance(segment, list):
                segment = str(segment)[1:-1].replace(' ', '').replace("'", '')
            segment = f'{{{segment}}}'
        # Extract and format sample
        sample = kwargs.get('sample', '')
        if sample:
            sample = f'@{sample}s'
        # Populate primekeys dict with Time and Wavelength values
        if start_time and end_time:
            # Check whether any primekey listed in PRIMEKEY_LIST_TIME has been passed through
            # PrimeKey() attribute. If yes, raise an error, since Time can only be passed
            # either through PrimeKey() attribute or Time() attribute.
            if all(x not in PRIMEKEY_LIST_TIME for x in primekey):
                timestr = '{start}-{end}{sample}'.format(
                    start=start_time.tai.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                    end=end_time.tai.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                    sample=sample)
            else:
                error_message = ("Time attribute has been passed both as a Time()"
                                 " and PrimeKey(). Please provide any one of them"
                                 " or separate them by OR operator."
                                 )
                raise ValueError(error_message)
        else:
            # This is executed when Time has not been passed through Time() attribute.
            # `match` stores all the time-type prime-keys that has been passed through
            # PrimeKey() attribute. The length of `match` won't ever be greater than 1,
            # but it is a good idea to keep a check.
            match = set(primekey.keys()) & PRIMEKEY_LIST_TIME
            if len(match) > 1:
                error_message = ("Querying of series, having more than 1 Time-type "
                                 "prime-keys is not yet supported. Alternative is to "
                                 "use only one of the primekey to query for series data.")
                raise ValueError(error_message)
            timestr = f"{primekey.pop(list(match)[0], '')}" if match else ''
        if wavelength == '':
            # This is executed when wavelength has been passed through PrimeKey().
            wavelength = f"{primekey.pop('WAVELNTH', '')}"
        elif not primekey.get('WAVELNTH'):
            wavelength = f'{int(np.ceil(wavelength.to_value(u.AA)))}'
        else:
            # This is executed when wavelength has been passed both through PrimeKey()
            # and Wavelength().
            error_message = ("Wavelength attribute has been passed both as a Wavelength()"
                             " and PrimeKey(). Please provide any one of them"
                             " or separate them by OR operator.")
            raise ValueError(error_message)
        # Populate primekey dict with formatted Time and Wavlength.
        if timestr:
            primekey['TIME'] = timestr
        if wavelength != '':
            primekey['WAVELNTH'] = wavelength
        # Extract and format primekeys
        primekey_string = ''
        c = drms.Client()
        si = c.info(series)
        primekeys_istime = si.keywords.loc[si.primekeys].is_time
        for pkey in primekeys_istime.index.values:
            if len(primekey) <= 0:
                break
            if primekeys_istime[pkey]:
                primekey_string += f"[{primekey.pop('TIME', '')}]"
            else:
                primekey_string += f"[{primekey.pop(pkey, '')}]"
        if not primekey_string:
            # primekey_string cannot be totally empty
            #
            # Note that whilst it is technically possible to just search by series,
            # this is not allowed here, because some of these would be very large
            # searches that would make JSOC unhappy with sunpy
            raise ValueError("Time, Wavelength or an explicit PrimeKey must be specified.")
        keys = []
        keyword_info = c.keys(series)
        for key, value in keyword.items():
            if key not in keyword_info:
                raise ValueError(f"Keyword: '{key}' is not supported by series: {series}")
            # We have to ensure that any input values that are strings are quoted
            keys.append(f"{key}{value['operator']}'{value['value']}'")
        keyword_string = f"[? {' AND '.join(keys)} ?]" if keys else ""
        return f"{series}{primekey_string}{keyword_string}{segment}"

    def _lookup_records(self, iargs):
        """
        Do a LookData request to JSOC to workout what results the query returns.
        """
        is_meta = iargs.get('meta', False)
        keywords = iargs.get('keys', '**ALL**')
        client = drms.Client()
        if 'series' not in iargs:
            error_message = "Series must be specified for a JSOC Query"
            raise ValueError(error_message)
        if not isinstance(keywords, list) and not isinstance(keywords, str):
            error_message = "Keywords can only be passed as a list or comma-separated strings."
            raise TypeError(error_message)
        # Get a set of the PrimeKeys that exist for the given series, and check
        # whether the passed PrimeKeys is a subset of that.
        primekeys = client.pkeys(iargs['series'])
        primekeys_passed = iargs.get('primekey', None)  # primekeys_passes is a dict, with key-value pairs.
        if primekeys_passed is not None and not set(list(primekeys_passed.keys())) <= set(primekeys):
            error_message = f"Unexpected PrimeKeys were passed. The series {iargs['series']} supports the following Keywords: {primekeys}"
            raise ValueError(error_message.format(series=iargs['series'], primekeys=primekeys))
        # Raise special error for wavelength (even though the code would ignore it anyway)
        wavelength = iargs.get('wavelength', '')
        if wavelength != '' and 'WAVELNTH' not in primekeys:
            error_message = (f"The series {iargs['series']} does not support wavelength attribute. "
                             f"The following primekeys are supported {primekeys}")
            raise TypeError(error_message.format(series=iargs['series'], pkeys=primekeys))
        # Raise errors for segments
        # Get a set of the segments that exist for the given series, and check
        # whether the passed segments is a subset of that.
        series = client.info(iargs['series'])
        # Fetches all valid segment names
        segments = list(series.segments.index.values)
        segments_passed = iargs.get('segment', None)
        if segments_passed is not None:
            if not isinstance(segments_passed, list) and not isinstance(segments_passed, str):
                error_message = "Segments can only be passed as a comma-separated string or a list of strings."
                raise TypeError(error_message)
            elif isinstance(segments_passed, str):
                segments_passed = segments_passed.replace(' ', '').split(',')
            if not set(segments_passed) <= set(segments):
                error_message = f"Unexpected Segments were passed. The series {iargs['series']} contains the following Segments {segments}"
                raise ValueError(error_message.format(series=iargs['series'], segs=segments))
            iargs['segment'] = segments_passed
        # If Time has been passed as a PrimeKey, convert the Time object into TAI time scale,
        # and then, convert it to datetime object.
        ds = self._make_recordset(**iargs)
        # Convert the list of keywords into comma-separated string.
        if isinstance(keywords, list):
            key = str(keywords)[1:-1].replace(' ', '').replace("'", '')
        else:
            key = keywords
        log.debug(f"Running following query: {ds}")
        log.debug(f"Requesting following keywords: {key}")
        result = client.query(ds, key=key, seg=segments_passed, rec_index=is_meta)
        if segments_passed:
            # DRMS returns a tuple of length 2 if there is segment information
            query_result, segment_result = result[0], result[1]
            query_result = pd.concat([query_result, segment_result],axis=1)
        else:
            query_result = result
        if result is None or query_result.empty:
            return astropy.table.Table()
        return astropy.table.Table.from_pandas(query_result)

    @classmethod
    def _can_handle_query(cls, *query):
        # Import here to prevent circular imports
        from sunpy.net import attrs as a

        required = {a.jsoc.Series}
        optional = {
            a.jsoc.Protocol, a.jsoc.Notify, a.Wavelength, a.Time,
            a.jsoc.Segment, a.jsoc.PrimeKey, a.Sample,
            a.jsoc.Cutout, a.jsoc.Keyword, a.jsoc.KeywordComparison,
        }
        return cls.check_attr_types_in_query(query, required, optional)

    @classmethod
    def _attrs_module(cls):
        return 'jsoc', 'sunpy.net.jsoc.attrs'

    @classmethod
    def register_values(cls):
        # We always use the local file for now.
        return cls.load_jsoc_values()

    @staticmethod
    def create_parse_jsoc_values():
        """
        Makes a network call to the VSO API that returns what keywords they support.
        We take this list and register all the keywords as corresponding Attrs.
        """
        here = os.path.dirname(os.path.realpath(__file__))
        client = drms.Client()
        # Series we are after
        data_sources = ["hmi", "mdi", "aia"]
        # Now get all the information we want.
        series_store = []
        segments = []
        for series in data_sources:
            info = client.series(rf'{series}\.')
            for item in info:
                try:
                    print(f'🛰 Getting info for {series}: {item}')
                    data = client.info(item)
                    series_store.append((data.name, data.note))
                    if not data.segments.empty:
                        segments.extend((row[0], row[1].iloc[-1]) for row in data.segments.iterrows())
                except Exception as e:
                    print(f"⚠️  {series} failed with error: {e}")
                    if item in ["hmi.V_avg120", "mdi.fdV_avg120"]:
                        # The following is from a private email from JSOC:
                        #
                        # So basically the 0xC5 LATIN-1 byte (Angstrom) got copied over, as-is, to the
                        # hmidb2 database. So hmidb2 is basically corrupt in a strict sense, although in
                        # practice, the only thing that gets messed up is this series, sometimes,
                        # because it is both slony-replicated and hmidb has a non-ascii char in it.
                        #
                        # Then jsocexintinfo.py asks PostgreSQL to "give me info for hmi.v_avg120, first
                        # converting to utf8 from ascii". PostgreSQL ignores any byte that is not an
                        # ascii byte and simply returns it as is. So now jsocexintfo.py has a byte,
                        # 0xC5, that it expects is a UTF-8 byte, which it isn't.
                        print(f"🛈 {item} has a known issue with the JSOC database.")
        series_store = list(set(series_store))
        segments = list(set(segments))
        with open(os.path.join(here, 'data', 'attrs.json'), 'w') as attrs_file:
            keyword_info = {
                "series_store": sorted(series_store),
                "segments": sorted(segments),
            }
            json.dump(keyword_info, attrs_file, indent=2)

    @staticmethod
    def load_jsoc_values():
        """
        We take this list and register all the keywords as corresponding Attrs.

        Returns
        -------
        dict
            The constructed Attrs dictionary ready to be passed into Attr registry.
        """
        # Import here to prevent circular imports
        from sunpy.net import attrs as a

        here = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(here, 'data', 'attrs.json')) as attrs_file:
            keyword_info = json.load(attrs_file)
        # Create attrs out of them.
        series_dict = {a.jsoc.Series: keyword_info["series_store"]}
        segments_dict = {a.jsoc.Segment: keyword_info["segments"]}
        attrs = series_dict | segments_dict
        return attrs
