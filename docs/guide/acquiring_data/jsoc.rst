---------------------------------------
Querying and Downloading Data from JSOC
---------------------------------------

Joint Science Operations Center (JSOC) contains data products from the Solar Dynamics Observatory,
as well as certain other missions and instruments. These data are available from the JSOC database,
which can be directly accessed by the online JSOC interface <http://jsoc.stanford.edu/ajax/lookdata.html>

SunPy's JSOC Client provides an easier interface to query for JSOC data and make export requests.
It uses drms module <https://github.com/kbg/drms> as its backend, and exposes a similar API as
the VSO Client.

Setup
-----

SunPy's JSOC Client is in ``sunpy.net``.  It can be imported as follows:

    >>> from sunpy.net import jsoc
    >>> client = jsoc.JSOCClient()

This creates your client object. The JSOC stages data before you can download it,
so a JSOC query is a three stage process, first you query the JSOC for records,
a table of these records is returned. Then you can request these records to be
staged for download and then you can download them.
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
Wavelength is another pre-defined attributes which is a prime-key.
Other primekeys which needs to be passed, should be manually passed in PrimeKey(). This
will be explained later in detail.

Constructing a Basic Query
^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's start with a very simple query.  We could ask for all `hmi.v_45s` series data
between January 1st and 2nd, 2014.

    >>> res = client.query(jsoc.attrs.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), jsoc.attrs.Series('hmi.v_45s'))

The variable ``qr`` is a Python list of
response objects, each one of which is a record found by the VSO. You can find how many
records were found by typing

    >>> len(res)
    81

To get a little bit more information about the records found, try

    >>> print(qr) # doctest:+SKIP
    ...


Now, let's break down the arguments of ``client.query`` to understand
better what we've done.  The first argument:

    ``jsoc.attrs.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00')``

sets the start and end times for the query (any date/time
format understood by SunPy's :ref:`parse_time function <parse-time>`
can be used to specify dates and time). The Time attribute takes UTC time,
as default. If you need to pass a Time in some other time scale, such as TAI,
pass an Astropy Time object, like:

	>>> from astropy.time import Time as T
	``jsoc.attrs.Time(T('2014-01-01T00:00:00', scale='tai'), T('2014-01-01T01:00:00', scale='tai'))``

The second argument:

    ``jsoc.attrs.Series('hmi.v_45s')``

sets the series we are looking for.

So what is going on here?
The notion is that a JSOC query has a set of attribute objects -
described in ``jsoc.attrs`` - that are specified to construct the query.

``jsoc.attrs.Series()`` is compulsory to be provided in each of the jsoc queries. Apart from this,
atleast one prime-key must be passed (generally ``jsoc.attrs.Time()``).

Querying with other PrimeKeys
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Other than Time, one other PrimeKey is supported with in-built attribute.
In case of AIA series, ``jsoc.attrs.Wavelength()`` can be passed as a prime-key.

    >>> import astropy.units as u	
	>>> res = client.query(jsoc.attrs.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
						   jsoc.attrs.Series('aia.lev1_euv_12s'), jsoc.attrs.Wavelength(304*u.AA))

Note that, only Time and Wavelength are in-built attributes here. If you need to pass any other primekey,
it should be passed like this:

	``jsoc.attrs.PrimeKey('HARPNUM', '4864')``

	or, if 2 or more PrimeKeys need to be passed together:
	``jsoc.attrs.PrimeKey('HARPNUM', '4864') & jsoc.attrs.PrimeKey('CAMERA', '2')``

Also, note that the pre-defined primkeys, Time and Wavelength can also be passed as above, but you need to
specify the exact keyword for it. For e.g. by :

	``jsoc.attrs.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), jsoc.attrs.PrimeKey('WAVELNTH', '161')``

If the correct keyword is not specified, or the passed PrimeKey is not supported by the given series, a
meaningful error will be thrown, which will give you the primekeys supported by that series. Hence, by looking
at the error, one can easily retry building the query with correct PrimeKeys.

Other important thing to note is that, Wavelength when passed through in-built attribute, should be passed as a
astropy quantity. Specifying spectral units in arguments is necessary or an error will be raised.
To know more check `astropy.units`.
But, when the same is passed through PrimeKey attribute, it should be passed as a string. All
other PrimeKey values passed through PrimeKey attribute, must be passed as a string.


Manually specifying keyword data to fetch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Upon doing ``client.query()`` as described above, only a limited set of keywords are returned in the response
object. These default keywords are 'DATE', 'TELESCOP', 'INSTRUME', 'T_OBS' and 'WAVELNTH'.

If you want to get a manual set of keywords in the response object, you can pass the set of keywords using
``jsoc.attrs.Keys()`` attribute.

	>>> res = client.query(jsoc.attrs.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
					       jsoc.attrs.Series('hmi.v_45s'),
					       jsoc.attrs.Keys(['TELESCOP', 'INSTRUME', 'T_OBS']))

The parameter passed into ``jsoc.attrs.Keys()`` can be either a list of strings, or a string with keywords seperated by
comma and a space. Meaning to say,

	``jsoc.attrs.Keys(['TELESCOP', 'INSTRUME', 'T_OBS'])`` and ``jsoc.attrs.Keys('TELESCOP, INSTRUME, T_OBS')``

both are correct.

Passing an incorrect keyword won't through an error, but the corresponding column in the astropy table will
contain ``Invalid KeyLink``.

To get all of the keywords, you can either use the ``search_metadata()`` method, or alternatively pass
``jsoc.attrs.Keys('***ALL***')`` with the series name and prime-key.


Using Segments
^^^^^^^^^^^^^^
In some cases, more than 1 file are present for the same set of query. These data are distinguished by what are called
`Segments`. It is necessary to specify the "Segment" which you need to download. Providing a segment won't have any affect
on the response object returned, but this will be required later, while making an export request.

A list of supported segments of a series, say ``hmi.sharp_720s`` can be obtained by :

	>>> import drms
	>>> c = drms.Client()
	>>> si = c.info('hmi.sharp_720s')
	>>> print(si.segments.index.values)

Also, if you provide an incorrect segment name, it will throw a meaningful error, specifying which segment values are supported
by the given series.

	>>> response = client.query(jsoc.attrs.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
                                jsoc.attrs.Series('aia.lev1_euv_12s'),
                                jsoc.attrs.Segment('image'))

To get files for more than 1 segment at the same time, chain ``jsoc.attrs.Segment()`` using ``AND`` operator.

	>>> res = client.query(jsoc.attrs.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
						   jsoc.attrs.Series('hmi.sharp_720s'),
						   jsoc.attrs.Segment('continuum') & jsoc.attrs.Segment('magnetogram'))


Using Sample
^^^^^^^^^^^^


Constructing complex queries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Complex queries can be built using OR operators.

Let's look for 2 dfifferent series data at the same time:

    >>> res = client.query(jsoc.attrs.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    					   jsoc.attrs.Series('hmi.v_45s') | jsoc.attrs.Series('aia.lev1_euv_12s'))
    >>> len(res)
    2188

The two series names are joined together by the operator "|".
This is the ``OR`` operator.  Think of the above query as setting a set
of conditions which get passed to the JSOC.  Let's say you want all the
EIT data from two separate days:

    >>> res = client.query(jsoc.attrs.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00') | 
                           jsoc.attrs.Time('2014-01-02T00:00:00', '2014-01-02T01:00:00'),
                           jsoc.attrs.Series('hmi.v_45s'))

Each of the arguments in this query style can be thought of as
setting conditions that the returned records must satisfy.

It should be noted that ``AND`` operator is supported by some of the attributes only. The attributes which
support "&" are ``PrimeKey()``, ``Segment()``. Using "&" with any other attributes will throw an error.