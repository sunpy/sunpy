***************************************
Querying and Downloading Data from JSOC
***************************************

Joint Science Operations Center (JSOC) contains data products from the Solar Dynamics Observatory,
as well as certain other missions and instruments. These data are available from the JSOC database,
which can be directly accessed by the online `JSOC interface <http://jsoc.stanford.edu/ajax/lookdata.html>`_.

SunPy's JSOC Client provides an easier interface to query for JSOC data and make export requests.
It uses `drms module <https://docs.sunpy.org/projects/drms>`_ as its backend, and exposes a similar API as
the VSO Client.

There are two ways of downloading JSOC data. One way is using Sunpy's unified search interface,
known as Fido. Fido supplies a single, easy and consistent way to to obtain most forms of solar physics data.
An alternative way to fetch data from JSOC is by using the underlying JSOC Client. This option
can be preferred when the complex searches are to be made, or when you need to separate the staging
and downloading steps, which is not supported by Fido.

The JSOC stages data before you can download it,
so a JSOC query is a three stage process. First you query the JSOC for records and
a table of these records is returned. Then you can request these records to be
staged for download and then you can download them. Fido combines the stages into 2,
`~sunpy.net.fido_factory.UnifiedDownloaderFactory.search` and
`~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`.

Setup
*****

SunPy's Fido module is in `sunpy.net`. It can be imported as follows:

    >>> from sunpy.net import Fido, attrs as a

The JSOC client handles the particulars of how the data from
the data provider is downloaded to your computer.

.. warning::

    You must have an email address registered with JSOC before you are allowed to make a request.
    See `this <http://jsoc.stanford.edu/ajax/register_email.html>`__ to register your email address.

Querying the JSOC
*****************

To search for data in JSOC, your query needs at minimum, a Series name and a PrimeKey.
Different PrimeKeys are supported by different Series, and you can find out the PrimeKeys
supported in any Series by:

    >>> import drms
    >>> c = drms.Client()  # doctest: +REMOTE_DATA
    >>> print(c.pkeys('hmi.m_720s'))  # doctest: +REMOTE_DATA
    ['T_REC', 'CAMERA']

The most common PrimeKey, that is supported by every Series is Time, that is denoted by
T_REC or T_OBS. Hence, Time can always be passed as an attribute while building a query.
Wavelength is another pre-defined attribute which is a PrimeKey.
Other PrimeKeys which need to be passed should be manually passed in
`~sunpy.net.jsoc.attrs.PrimeKey`. This will be explained later in detail.

Constructing a Basic Query
==========================

Let's start with a very simple query.  We could ask for all ``hmi.v_45s`` series data
between January 1st from 00:00 to 01:00, 2014.
We can add email address to the search query with the `sunpy.net.jsoc.attrs.Notify` attribute.
Please note you can search without this but right now, you can not add the email address after the search.

    >>> res = Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Series('hmi.v_45s'),
    ...                   a.jsoc.Notify('sunpy@sunpy.org'))  # doctest: +REMOTE_DATA

This returns an `~sunpy.net.fido_factory.UnifiedResponse` object containing
information on the available online files which fit the criteria specified by
the attrs objects in the above call. It does not download the files.

To see a summary of results of our query, simply type the name of the
variable set to the Fido search, in this case, res::

    >>> res  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    81 Results from the JSOCClient:
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
                        ...      ...        ...      ...     ...
    2014.01.01_00:54:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:54:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:55:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:56:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:58:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    Length = 81 rows
    <BLANKLINE>
    <BLANKLINE>


Now, let's break down the arguments of ``Fido.search`` to understand
better what we've done.  The first argument ``a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00')``
sets the start and end times for the query (any date/time
format understood by SunPy's :ref:`parse_time function <parse-time>`
can be used to specify dates and time). The Time attribute takes UTC time,
as default. If you need to pass a Time in some other time scale, such as TAI,
pass an Astropy Time object, like::

    >>> import astropy.time

Then, the Time attribute can be passed as::

    >>> a.Time(astropy.time.Time('2014-01-01T00:00:00', scale='tai'), astropy.time.Time('2014-01-01T01:00:00', scale='tai'))
    <sunpy.net.attrs.Time(2014-01-01 00:00:00.000, 2014-01-01 01:00:00.000)>

The second argument::

    >>> a.jsoc.Series('hmi.v_45s')
    <sunpy.net.jsoc.attrs.Series(hmi.v_45s: Dopplergrams with a cadence of 45 seconds) object ...>

sets the series we are looking for.

So what is going on here?
The notion is that a JSOC query has a set of attribute objects, imported as ``a.jsoc``,
that are specified to construct the query.

``a.jsoc.Series()`` is compulsory to be provided in each of the jsoc queries. Apart from this,
at least one PrimeKey must be passed (generally ``a.Time()``).

The third argument::

    >>> a.jsoc.Notify('sunpy@sunpy.org')
    <sunpy.net.jsoc.attrs.Notify: sunpy@sunpy.org object ...>

tells JSOC what email address you are registered with and to email when your request is ready to download.

Querying with other PrimeKeys
=============================

Other than Time, one other PrimeKey is supported with in-built attribute.
In case of AIA series, ``a.Wavelength()`` can be passed as a PrimeKey::

    >>> import astropy.units as u
    >>> res = Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                               a.jsoc.Notify('sunpy@sunpy.org'),
    ...                               a.jsoc.Series('aia.lev1_euv_12s'),
    ...                               a.Wavelength(304*u.AA))  # doctest: +REMOTE_DATA

Note that, only Time and Wavelength are in-built attributes here. If you need to pass any other PrimeKey,
it should be passed like this::

    >>> a.jsoc.PrimeKey('HARPNUM', '4864')
    <sunpy.net.jsoc.attrs.PrimeKey object at ...>
    ('HARPNUM', '4864')

If 2 or more PrimeKeys need to be passed together::

    >>> a.jsoc.PrimeKey('HARPNUM', '4864') & a.jsoc.PrimeKey('CAMERA', '2')
    <AttrAnd([<sunpy.net.jsoc.attrs.PrimeKey object at ...>
    ('HARPNUM', '4864'), <sunpy.net.jsoc.attrs.PrimeKey object at ...>
    ('CAMERA', '2')])>

Also, note that the pre-defined primekeys, Time and Wavelength can also be passed as above, but you need to
specify the exact keyword for it. For e.g. by::

    >>> a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.PrimeKey('WAVELNTH', '161')
    (<sunpy.net.attrs.Time(2014-01-01 00:00:00.000, 2014-01-01 01:00:00.000)>, <sunpy.net.jsoc.attrs.PrimeKey object at ...>
    ('WAVELNTH', '161'))

If the correct keyword is not specified, or the passed PrimeKey is not supported by the given series, a
meaningful error will be thrown, which will give you the PrimeKeys supported by that series. Hence, by looking
at the error, one can easily retry building the query with correct PrimeKeys.

Another important thing to note is that, Wavelength when passed through in-built attribute, should be passed as an
Astropy quantity. Specifying spectral units in arguments is necessary or an error will be raised.
For more information on units, see `~astropy.units`.
But, when the same is passed through PrimeKey attribute, it should be passed as a string. All
other PrimeKey values passed through PrimeKey attribute, must be passed as a string.


Manually specifying keyword data to fetch
=========================================

Upon doing ``Fido.search()`` as described above, only a limited set of keywords are returned in the response
object. These default keywords are ``'DATE'``, ``'TELESCOP'``, ``'INSTRUME'``, ``'T_OBS'`` and ``'WAVELNTH'``.

If you want to get a manual set of keywords in the response object, you can pass the set of keywords using
:meth:`~sunpy.net.base_client.BaseQueryResponseTable.show` method.

    >>> res = Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                   a.jsoc.Series('hmi.v_45s'), a.jsoc.Notify('sunpy@sunpy.org'))  # doctest: +REMOTE_DATA
    >>> res.show('TELESCOP', 'INSTRUME', 'T_OBS')  # doctest: +REMOTE_DATA
    [<QTable length=81>
    TELESCOP  INSTRUME           T_OBS
    str7     str10             str23
    -------- ---------- -----------------------
    SDO/HMI HMI_FRONT2 2014.01.01_00:00:37_TAI
    SDO/HMI HMI_FRONT2 2014.01.01_00:01:22_TAI
    SDO/HMI HMI_FRONT2 2014.01.01_00:02:07_TAI
    SDO/HMI HMI_FRONT2 2014.01.01_00:02:52_TAI
        ...        ...                     ...
    SDO/HMI HMI_FRONT2 2014.01.01_00:57:37_TAI
    SDO/HMI HMI_FRONT2 2014.01.01_00:58:22_TAI
    SDO/HMI HMI_FRONT2 2014.01.01_00:59:07_TAI
    SDO/HMI HMI_FRONT2 2014.01.01_00:59:52_TAI
    SDO/HMI HMI_FRONT2 2014.01.01_01:00:37_TAI]

Passing an incorrect keyword won't throw an error, but the corresponding column in the table will
not be displayed.

To display all of the columns, we can use ``show()`` without passing any arguments::

    >>> res.show()  # doctest: +REMOTE_DATA
    [<QTable length=81>
            DATE                DATE__OBS        ... CALVER64
           str20                  str23          ...  int64
    -------------------- ----------------------- ... --------
    2014-01-05T17:46:02Z 2013-12-31T23:59:39.20Z ...     4370
    2014-01-05T17:47:10Z 2014-01-01T00:00:24.20Z ...     4370
    2014-01-05T17:48:18Z 2014-01-01T00:01:09.20Z ...     4370
    2014-01-05T17:49:25Z 2014-01-01T00:01:54.20Z ...     4370
                     ...                     ... ...      ...
    2014-01-05T17:41:25Z 2014-01-01T00:56:39.20Z ...     4370
    2014-01-05T17:42:33Z 2014-01-01T00:57:24.20Z ...     4370
    2014-01-05T17:43:41Z 2014-01-01T00:58:09.20Z ...     4370
    2014-01-05T17:44:52Z 2014-01-01T00:58:54.20Z ...     4370
    2014-01-05T17:46:03Z 2014-01-01T00:59:39.20Z ...     4370]

or you can print the results::

    >>> res  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    81 Results from the JSOCClient:
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
                        ...      ...        ...      ...     ...
    2014.01.01_00:53:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:54:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:54:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:55:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:56:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:58:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    Length = 81 rows
    <BLANKLINE>
    <BLANKLINE>

Using Segments
==============
In some cases, more than 1 file are present for the same set of query. These data are distinguished by what are called
Segments. It is necessary to specify the "Segment" which you need to download. Providing a segment won't have any affect
on the response object returned, but this will be required later, while making an export request.

A list of supported segments of a series, say ``hmi.sharp_720s`` can be obtained by::

    >>> import drms
    >>> c = drms.Client()  # doctest: +REMOTE_DATA
    >>> si = c.info('hmi.sharp_720s')  # doctest: +REMOTE_DATA
    >>> print(si.segments.index.values)  # doctest: +REMOTE_DATA
    ['magnetogram' 'bitmap' 'Dopplergram' 'continuum' 'inclination' 'azimuth'
     'field' 'vlos_mag' 'dop_width' 'eta_0' 'damping' 'src_continuum'
     'src_grad' 'alpha_mag' 'chisq' 'conv_flag' 'info_map' 'confid_map'
     'inclination_err' 'azimuth_err' 'field_err' 'vlos_err' 'alpha_err'
     'field_inclination_err' 'field_az_err' 'inclin_azimuth_err'
     'field_alpha_err' 'inclination_alpha_err' 'azimuth_alpha_err' 'disambig'
     'conf_disambig']

Also, if you provide an incorrect segment name, it will throw a meaningful error, specifying which segment values are supported
by the given series::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...             a.jsoc.Series('hmi.sharp_720s'), a.jsoc.Notify('sunpy@sunpy.org'),
    ...             a.jsoc.Segment('image'))  # doctest: +REMOTE_DATA
    Traceback (most recent call last):
    ...
    ValueError: Unexpected Segments were passed. The series hmi.sharp_720s contains the following Segments ['magnetogram', 'bitmap', 'Dopplergram', 'continuum', 'inclination', 'azimuth', 'field', 'vlos_mag', 'dop_width', 'eta_0', 'damping', 'src_continuum', 'src_grad', 'alpha_mag', 'chisq', 'conv_flag', 'info_map', 'confid_map', 'inclination_err', 'azimuth_err', 'field_err', 'vlos_err', 'alpha_err', 'field_inclination_err', 'field_az_err', 'inclin_azimuth_err', 'field_alpha_err', 'inclination_alpha_err', 'azimuth_alpha_err', 'disambig', 'conf_disambig']


To get files for more than 1 segment at the same time, chain ``a.jsoc.Segment()`` using ``AND`` operator::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...             a.jsoc.Series('hmi.sharp_720s'), a.jsoc.Notify('sunpy@sunpy.org'),
    ...             a.jsoc.Segment('continuum') & a.jsoc.Segment('magnetogram'))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    61 Results from the JSOCClient:
             T_REC          TELESCOP  INSTRUME WAVELNTH CAR_ROT
    ----------------------- -------- --------- -------- -------
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
                        ...      ...       ...      ...     ...
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:12:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:24:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:36:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_00:48:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2145
    Length = 61 rows
    <BLANKLINE>
    <BLANKLINE>


Using Sample
============
In case you need to query for data, at some interval of time, say every 10 min, you can pass it
using `~sunpy.net.attrs.Sample`. In other words, if you need to query for ``hmi.v_45s`` series data
between January 1st from 00:00 to 01:00, 2014, every 10 minutes, you can do::

    >>> import astropy.units as u
    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Notify('sunpy@sunpy.org'),
    ...             a.jsoc.Series('hmi.v_45s'), a.Sample(10*u.min))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    7 Results from the JSOCClient:
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
    ----------------------- -------- ---------- -------- -------
    2014.01.01_00:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:10:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:20:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:30:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:39:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:49:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    <BLANKLINE>
    <BLANKLINE>

Note that the argument passed in ``a.Sample()`` must be an Astropy quantity, convertible
into seconds.

Constructing complex queries
============================

Complex queries can be built using ``OR`` operators.

Let's look for 2 different series data at the same time::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Notify('sunpy@sunpy.org'),
    ...             a.jsoc.Series('hmi.v_45s') | a.jsoc.Series('aia.lev1_euv_12s'))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    81 Results from the JSOCClient:
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
                        ...      ...        ...      ...     ...
    2014.01.01_00:54:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:54:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:55:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:56:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:58:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    Length = 81 rows
    <BLANKLINE>
    2107 Results from the JSOCClient:
           T_REC         TELESCOP INSTRUME WAVELNTH CAR_ROT
    -------------------- -------- -------- -------- -------
    2014-01-01T00:00:01Z  SDO/AIA    AIA_4       94    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_1      131    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_2      193    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_2      211    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T00:00:01Z  SDO/AIA    AIA_1      335    2145
    2014-01-01T00:00:13Z  SDO/AIA    AIA_4       94    2145
    2014-01-01T00:00:13Z  SDO/AIA    AIA_1      131    2145
    2014-01-01T00:00:13Z  SDO/AIA    AIA_3      171    2145
                     ...      ...      ...      ...     ...
    2014-01-01T00:59:49Z  SDO/AIA    AIA_2      211    2145
    2014-01-01T00:59:49Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T00:59:49Z  SDO/AIA    AIA_1      335    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_4       94    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_1      131    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_3      171    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_2      193    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_2      211    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_4      304    2145
    2014-01-01T01:00:01Z  SDO/AIA    AIA_1      335    2145
    Length = 2107 rows
    <BLANKLINE>
    <BLANKLINE>

The two series names are joined together by the operator ``|``.
This is the ``OR`` operator.  Think of the above query as setting a set
of conditions which get passed to the JSOC.  Let's say you want all the
``hmi.v_45s`` data from two separate days::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00') |
    ...             a.Time('2014-01-02T00:00:00', '2014-01-02T01:00:00'),
    ...             a.jsoc.Series('hmi.v_45s'), a.jsoc.Notify('sunpy@sunpy.org'))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    81 Results from the JSOCClient:
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
                        ...      ...        ...      ...     ...
    2014.01.01_00:54:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:54:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:55:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:56:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:57:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:58:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.01_01:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    Length = 81 rows
    <BLANKLINE>
    81 Results from the JSOCClient:
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
    ----------------------- -------- ---------- -------- -------
    2014.01.02_00:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:01:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:02:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:03:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:03:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:04:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:05:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:06:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:06:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:07:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
                        ...      ...        ...      ...     ...
    2014.01.02_00:54:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:54:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:55:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:56:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:57:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:57:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:58:30_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_00:59:15_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_01:00:00_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    2014.01.02_01:00:45_TAI  SDO/HMI HMI_FRONT2   6173.0    2145
    Length = 81 rows
    <BLANKLINE>
    <BLANKLINE>

Each of the arguments in this query style can be thought of as
setting conditions that the returned records must satisfy.

It should be noted that ``AND`` operator is supported by some of the attributes only. The attributes which
support "&" are `~sunpy.net.jsoc.attrs.PrimeKey` and `~sunpy.net.jsoc.attrs.Segment`.
Using "&" with any other attributes will throw an error.

Downloading data
****************

To download the files located by `~sunpy.net.fido_factory.UnifiedDownloaderFactory.search`,
you can download them by `~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`::

    >>> downloaded_files = Fido.fetch(res)  # doctest: +SKIP

Using JSOCClient for complex usage
**********************************

Fido interface uses `~sunpy.net.jsoc.JSOCClient` in its backend, and combines
the last 2 stages the JSOC process into one. You can directly use the JSOC
client to make queries, instead of the Fido client. This will allow you to
separate the 3 stages of the JSOC process, and perform it individually, hence
allowing a greater control over the whole process.

Setup
=====

SunPy's JSOC module is in `~sunpy.net`.  It can be imported as follows::

    >>> from sunpy.net import jsoc
    >>> client = jsoc.JSOCClient()  # doctest: +REMOTE_DATA

This creates your client object.


Making a query
==============

Querying JSOC using the JSOC client is very similar to what we were doing with Fido.
As above, we have to make sure we have an email address registered with JSOC before you are allowed to make a request.
See `this <http://jsoc.stanford.edu/ajax/register_email.html>`__ to register your email address.
We can add an email address to the search query with the `sunpy.net.jsoc.attrs.Notify` attribute.
Please note you can search without this but right now, you can not add the email address after the search::

    >>> from sunpy.net import attrs as a
    >>> res = client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                     a.jsoc.Series('hmi.v_45s'),
    ...                     a.jsoc.Notify('sunpy@sunpy.org'))  # doctest: +REMOTE_DATA

Apart from the function name, everything is the same. You need to pass the same values in the
`~sunpy.net.jsoc.JSOCClient.search` as you did in `~sunpy.net.fido_factory.UnifiedDownloaderFactory.search`.
Complex queries can be built in a similar way, and all other things are the same.

Staging the request
===================

JSOC is a 3-stage process, and after getting the query results, we need to stage a request for the data to be
downloaded. Only then, can we download them. The download request can be staged like this::

    >>> requests = client.request_data(res)  # doctest: +SKIP
    >>> print(requests)  # doctest: +SKIP
    <ExportRequest id="JSOC_20170713_1461", status=0>

The function `~sunpy.net.jsoc.JSOCClient.request_data` stages the request.
It returns a `drms.ExportRequest` object, which has many attributes.
The most important ones are ``id`` and ``status``. Only when the status is 0, we can
move to the third step, i.e. downloading the data.

If you are making more than 1 query at a time, it will return a list of `~drms.ExportRequest` objects. Hence, access the
list elements accordingly. You can get the id and status of the request (if it is not a list) by::

    >>> requests.id  # doctest: +SKIP
    JSOC_20170713_1461
    >>> requests.status  # doctest: +SKIP
    0


Downloading data
================

Once the status code is 0 you can download the data using the
`~sunpy.net.jsoc.JSOCClient.get_request` method::

    >>> res = client.get_request(requests)  # doctest: +SKIP

This returns a Results instance which can be used to watch the progress of the download::

    >>> res.wait(progress=True)   # doctest: +SKIP
