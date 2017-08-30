---------------------------------------
Querying and Downloading Data from JSOC
---------------------------------------

Joint Science Operations Center (JSOC) contains data products from the Solar Dynamics Observatory,
as well as certain other missions and instruments. These data are available from the JSOC database,
which can be directly accessed by the online `JSOC interface <http://jsoc.stanford.edu/ajax/lookdata.html>`

SunPy's JSOC Client provides an easier interface to query for JSOC data and make export requests.
It uses `drms module <https://github.com/kbg/drms>` as its backend, and exposes a similar API as
the VSO Client.

There are two ways of downloading JSOC data. One way is using Sunpy's unified search interface,
known as Fido. Fido supplies a single, easy and consistent way to to obtain most forms of physics
data. An alternative way to fetch data from JSOC is by using the underlying JSOC Client. This option
can be preferred when the complex searches are to be made, or when you need to separate the staging
and downloading steps, which is not supported by Fido.

The JSOC stages data before you can download it,
so a JSOC query is a three stage process, first you query the JSOC for records,
a table of these records is returned. Then you can request these records to be
staged for download and then you can download them. Fido combines the stages into 2,
`~sunpy.net.fido_factory.UnifiedDownloaderFactory.search` and
`~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`.

Setup
-----

SunPy's FIDO module is in `sunpy.net`.  It can be imported as follows:

    >>> from sunpy.net import Fido, attrs as a

The JSOC client handles the particulars of how the data from
the data provider is downloaded to your computer.

Querying the JSOC
-----------------

To search for data in JSOC, your query needs at minimum a series name and a primekey.
Different primekeys are supported by different series, and you can find out the primekeys
supported in any series by:

	>>> import drms
	>>> c = drms.Client()
	>>> print(c.pkeys('hmi.m_720s'))

The most common primekey, that is supported by every series is Time, that is denoted by
T_REC or T_OBS. Hence, Time can always be passed as an attribute while building a query.
Wavelength is another pre-defined attributes which is a primekey.
Other primekeys which needs to be passed, should be manually passed in
`~sunpy.net.jsoc.attrs.PrimeKey`. This will be explained later in detail.

Constructing a Basic Query
^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's start with a very simple query.  We could ask for all ``hmi.v_45s`` series data
between January 1st from 00:00 to 01:00, 2014.

    >>> res = Fido.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Series('hmi.v_45s'))

This returns an `~sunpy.net.fido_factory.UnifiedResponse` object containing
information on the available online files which fit the criteria specified by
the attrs objects in the above call. It does not download the files.

To see a summary of results of our query, simple type the name of the
variable set to the Fido search, in this case, res::

    >>> res

Now, let's break down the arguments of ``Fido.search`` to understand
better what we've done.  The first argument::

    ``a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00')``

sets the start and end times for the query (any date/time
format understood by SunPy's :ref:`parse_time function <parse-time>`
can be used to specify dates and time). The Time attribute takes UTC time,
as default. If you need to pass a Time in some other time scale, such as TAI,
pass an Astropy Time object, like:

	>>> import astropy.time

Then, the Time attribute can be passed as::

	``a.jsoc.Time(astropy.time.Time('2014-01-01T00:00:00', scale='tai'),
	              astropy.time.Time('2014-01-01T01:00:00', scale='tai'))``

The second argument::

    ``a.jsoc.Series('hmi.v_45s')``

sets the series we are looking for.

So what is going on here?
The notion is that a JSOC query has a set of attribute objects -
described in ``a.jsoc`` - that are specified to construct the query.

``a.jsoc.Series()`` is compulsory to be provided in each of the jsoc queries. Apart from this,
atleast one primekey must be passed (generally ``a.jsoc.Time()``).

Querying with other PrimeKeys
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Other than Time, one other PrimeKey is supported with in-built attribute.
In case of AIA series, ``a.jsoc.Wavelength()`` can be passed as a primekey.

    >>> import astropy.units as u	
	>>> res = Fido.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
						  a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Wavelength(304*u.AA))

Note that, only Time and Wavelength are in-built attributes here. If you need to pass any other primekey,
it should be passed like this:

	``a.jsoc.PrimeKey('HARPNUM', '4864')``

	or, if 2 or more PrimeKeys need to be passed together:
	``a.jsoc.PrimeKey('HARPNUM', '4864') & a.jsoc.PrimeKey('CAMERA', '2')``

Also, note that the pre-defined primkeys, Time and Wavelength can also be passed as above, but you need to
specify the exact keyword for it. For e.g. by :

	``a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.PrimeKey('WAVELNTH', '161')``

If the correct keyword is not specified, or the passed PrimeKey is not supported by the given series, a
meaningful error will be thrown, which will give you the primekeys supported by that series. Hence, by looking
at the error, one can easily retry building the query with correct PrimeKeys.

Other important thing to note is that, Wavelength when passed through in-built attribute, should be passed as a
astropy quantity. Specifying spectral units in arguments is necessary or an error will be raised.
For more information on units, see `astropy.units`.
But, when the same is passed through PrimeKey attribute, it should be passed as a string. All
other PrimeKey values passed through PrimeKey attribute, must be passed as a string.


Manually specifying keyword data to fetch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Upon doing ``Fido.search()`` as described above, only a limited set of keywords are returned in the response
object. These default keywords are ``'DATE'``, ``'TELESCOP'``, ``'INSTRUME'``, ``'T_OBS'`` and ``'WAVELNTH'``.

If you want to get a manual set of keywords in the response object, you can pass the set of keywords using
`~sunpy.net.jsoc.attrs.Keys` attribute.

	>>> res = Fido.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
					      a.jsoc.Series('hmi.v_45s'),
					      a.jsoc.Keys(['TELESCOP', 'INSTRUME', 'T_OBS']))

The parameter passed into ``a.jsoc.Keys()`` can be either a list of strings, or a string with keywords seperated by
comma and a space. Meaning to say,::

	``a.jsoc.Keys(['TELESCOP', 'INSTRUME', 'T_OBS'])`` and ``jsoc.attrs.Keys('TELESCOP, INSTRUME, T_OBS')``

both are correct.

Passing an incorrect keyword won't through an error, but the corresponding column in the astropy table will
contain ``Invalid KeyLink``.

To get all of the keywords, you can either use the `~sunpy.net.jsoc.JSOCClient.search_metadata` method,
or alternatively pass ``a.jsoc.Keys('***ALL***')`` with the series name and primekey.


Using Segments
^^^^^^^^^^^^^^
In some cases, more than 1 file are present for the same set of query. These data are distinguished by what are called
``Segments``. It is necessary to specify the "Segment" which you need to download. Providing a segment won't have any affect
on the response object returned, but this will be required later, while making an export request.

A list of supported segments of a series, say ``hmi.sharp_720s`` can be obtained by :

	>>> import drms
	>>> c = drms.Client()
	>>> si = c.info('hmi.sharp_720s')
	>>> print(si.segments.index.values)

Also, if you provide an incorrect segment name, it will throw a meaningful error, specifying which segment values are supported
by the given series.

	>>> response = Fido.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
                               a.jsoc.Series('aia.lev1_euv_12s'),
                               a.jsoc.Segment('image'))

To get files for more than 1 segment at the same time, chain ``a.jsoc.Segment()`` using ``AND`` operator.

	>>> res = Fido.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
						  a.jsoc.Series('hmi.sharp_720s'),
						  a.jsoc.Segment('continuum') & a.jsoc.Segment('magnetogram'))


Using Sample
^^^^^^^^^^^^
In case you need to query for data, at some interval of time, say every 10 min, you can pass it
using `~sunpy.net.jsoc.attrs.Sample`. In other words, if you need to query for `hmi.v_45s` series data
between January 1st from 00:00 to 01:00, 2014, every 10 minutes, you can do:

	>>> import astropy.units as u
	>>> res = Fido.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
						  a.jsoc.Series('hmi.v_45s'), a.jsoc.Sample(10*u.min))

Note that the argument passed in ``a.jsoc.Sample()`` must be an astropy quanitity, convertible
into seconds.

Constructing complex queries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Complex queries can be built using OR operators.

Let's look for 2 dfifferent series data at the same time:

    >>> res = Fido.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    					  a.jsoc.Series('hmi.v_45s') | a.jsoc.Series('aia.lev1_euv_12s'))

The two series names are joined together by the operator "|".
This is the ``OR`` operator.  Think of the above query as setting a set
of conditions which get passed to the JSOC.  Let's say you want all the
hmi.v_45s data from two separate days:

    >>> res = Fido.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00') | 
                          a.jsoc.Time('2014-01-02T00:00:00', '2014-01-02T01:00:00'),
                          a.jsoc.Series('hmi.v_45s'))

Each of the arguments in this query style can be thought of as
setting conditions that the returned records must satisfy.

It should be noted that ``AND`` operator is supported by some of the attributes only. The attributes which
support "&" are `~sunpy.net.jsoc.attrs.PrimeKey` and `~sunpy.net.jsoc.attrs.Segment`.
Using "&" with any other attributes will throw an error.

Downloading data
----------------

To download the files located by `~sunpy.net.fido_factory.UnifiedDownloaderFactory.search`,
you can download them by `~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`:

	>>> downloaded_files = Fido.fetch(res)

Using JSOCClient for complex usage
----------------------------------

Fido interface uses JSOCClient in its backend, which combines the last 2 stages the JSOC process into
one. You can directly use the JSOC Client to make queries, instead of the Fido Client. This will allow you
to separate the 3 stages of the JSOC process, and perform it individually, hence allowing a greater
control over the whole process.

Setup
^^^^^

SunPy's JSOC module is in `sunpy.net`.  It can be imported as follows:

    >>> from sunpy.net import jsoc
    >>> client = jsoc.JSOCClient()

This creates your client object.

Making a query
^^^^^^^^^^^^^^

Querying JSOC using the JSOCClient is completely similar to what we were doing with Fido.

	>>> from sunpy.net import attrs as a
	>>> res = client.search(a.jsoc.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Series('hmi.v_45s'))

Apart from the function name, everything is same. You need to pass the same values in the
`~sunpy.net.jsoc.JSOCClient.search` as you did in `~sunpy.net.fido_factory.UnifiedDownloaderFactory.search`.
Complex queries can be built in a similar way, and all other things are same.

Staging the request
^^^^^^^^^^^^^^^^^^^

JSOC is a 3-stage process, and after getting the query results, we need to stage a request for the data to be
downloaded. Only then, can we download them. The download request can be staged like this:

	>>> requests = client.request_data(res)
	>>> print(requests)

	<ExportRequest id="JSOC_20170713_1461", status=0>

The function `~sunpy.net.jsoc.JSOCClient.request_data` stages the request.
It returns a `drms.ExportRequest` object, which has many attributes.
The most important ones are ``ExportRequest id`` and ``status``. Only when the status is 0, we can
move to the third step, i.e. downloading the data.

If you are making more than 1 query at a time, it will return a list of ExportRequest objects. Hence, access the
list elements accordingly. You can get the id and status of the request (if it is not a list) by:

	>>> requests.id
	>>> requests.status

You can also check the status of a request made by:

	>>> status = client.check_request(requests)

You can pass a list of ExportRequest objects, and a list of status' will be returned.

Downloading data
^^^^^^^^^^^^^^^^

Once the status code is 0 you can download the data using the
`~sunpy.net.jsoc.JSOCClient.get_request` method:

    >>> res = client.get_request(requests)

This returns a Results instance which can be used to watch the progress of the download.

    >>> res.wait(progress=True)   # doctest: +SKIP