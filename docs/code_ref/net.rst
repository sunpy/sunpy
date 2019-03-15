SunPy net
=========

SunPy's net submodule contains a lot of different code for accessing various
solar physics related web services. This submodule contains many layers. Most
users should use `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`, which
is an interface to multiple sources including all the sources implemented in
`~sunpy.net.dataretriever` as well as `~sunpy.net.vso` and `~sunpy.net.jsoc`.
Fido ~`sunpy.net.fido_factory.UnifiedDownloaderFactory` can be used like so::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Instrument('lyra'))  # doctest: +REMOTE_DATA
    >>> files = Fido.fetch(results)  # doctest: +SKIP

---------------------------------------
Finding and Downloading Data using Fido
---------------------------------------

This guide outlines how to search for and download data using SunPy's
Federated Internet Data Obtainer...or more usually (and sanely) referred to as Fido.
`Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` is a unified interface for searching
and fetching solar physics data irrespective of the underlying
client or webservice through which the data is obtained, e.g. VSO_,
JSOC_, etc.  It therefore supplies a single, easy and consistent way to
obtain most forms of solar physics data.

Import
------

The `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` object is in
`sunpy.net`. It can be imported as follows::

    >>> from sunpy.net import Fido, attrs as a

Searching for Data Using Fido
-----------------------------

To search for data with Fido, you need to specify attributes to search against.
The (partial) range of allowed attributes are found in the `vso.attrs <sunpy.net.vso.attrs>`
and `jsoc.attrs <sunpy.net.jsoc.attrs>`.
Examples of these attributes are `a.Time <sunpy.net.vso.attrs.Time>`,
`a.Instrument <sunpy.net.vso.attrs.Instrument>`,
`a.Wavelength <sunpy.net.vso.attrs.Wavelength>`, some of these attributes are
client specific, such as `a.vso <sunpy.net.vso.attrs>` or
`a.jsoc <sunpy.net.jsoc.attrs>`.::

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('lyra')) # doctest: +REMOTE_DATA

This returns an `~sunpy.net.fido_factory.UnifiedResponse` object containing
information on the available online files which fit the criteria specified by
the attrs objects in the above call. It does not download the files. For
instructions on how to download data using Fido, see :ref:`downloading_data`.

To see a summary of results of our query, simple type the name of the
variable set to the Fido search, in this case, result::

    >>> result  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the LYRAClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str6     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2012-03-04 00:00:00 2012-03-06 00:00:00 Proba2       lyra        nan
    2012-03-04 00:00:00 2012-03-06 00:00:00 Proba2       lyra        nan
    2012-03-04 00:00:00 2012-03-06 00:00:00 Proba2       lyra        nan
    <BLANKLINE>
    <BLANKLINE>

Queries can be made more flexible or specific by adding more attrs objects to
the `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` search. Specific
passbands can be searched for by supplying an `~astropy.units.Quantity` to the
`a.Wavelength <sunpy.net.vso.attrs.Wavelength>` attribute::

    >>> import astropy.units as u
    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('norh'),
    ...             a.Wavelength(17*u.GHz))  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the NoRHClient:
         Start Time           End Time      Source Instrument   Wavelength
           str19               str19         str4     str4        str14
    ------------------- ------------------- ------ ---------- --------------
    2012-03-04 00:00:00 2012-03-05 00:00:00   NAOJ       NORH 17000000.0 kHz
    2012-03-05 00:00:00 2012-03-06 00:00:00   NAOJ       NORH 17000000.0 kHz
    2012-03-06 00:00:00 2012-03-07 00:00:00   NAOJ       NORH 17000000.0 kHz
    <BLANKLINE>
    <BLANKLINE>

Data of a given cadence can also be specified using the Sample attribute. To
search for data at a given cadence use the
`a.vso.Sample <sunpy.net.vso.attrs.Sample>` attribute.
`a.vso.Sample <sunpy.net.vso.attrs.Sample>` is only supported by the
`sunpy.net.vso.VSOClient` hence it has the ``a.vso`` prefix. Attributes
like this which are client specific will result in
`Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` only searching that
client for results, in this case VSO.::

    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('aia'),
    ...             a.Wavelength(171*u.angstrom), a.vso.Sample(10*u.minute))  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    289 Results from the VSOClient:
       Start Time [1]       End Time [1]    Source ...   Type   Wavelength [2]
                                                   ...             Angstrom
           str19               str19         str3  ...   str8      float64
    ------------------- ------------------- ------ ... -------- --------------
    2012-03-05 04:30:00 2012-03-05 04:30:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-04 09:10:00 2012-03-04 09:10:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-05 22:50:00 2012-03-05 22:50:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-05 20:50:00 2012-03-05 20:50:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-04 01:00:00 2012-03-04 01:00:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-04 15:40:00 2012-03-04 15:40:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-05 12:40:00 2012-03-05 12:40:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-05 10:50:00 2012-03-05 10:50:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-04 01:50:00 2012-03-04 01:50:01    SDO ... FULLDISK 171.0 .. 171.0
                    ...                 ...    ... ...      ...            ...
    2012-03-04 07:40:00 2012-03-04 07:40:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-05 05:10:00 2012-03-05 05:10:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-05 08:00:00 2012-03-05 08:00:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-04 20:10:00 2012-03-04 20:10:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-04 07:30:00 2012-03-04 07:30:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-05 00:00:00 2012-03-05 00:00:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-05 13:50:00 2012-03-05 13:50:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-05 15:20:00 2012-03-05 15:20:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-05 08:10:00 2012-03-05 08:10:01    SDO ... FULLDISK 171.0 .. 171.0
    2012-03-04 02:50:00 2012-03-04 02:50:01    SDO ... FULLDISK 171.0 .. 171.0
    <BLANKLINE>
    <BLANKLINE>

To search for data from multiple instruments, wavelengths, times etc., use the
pipe ``|`` operator. This joins queries together just as the logical ``OR``
operator would::

    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'),
    ...             a.Instrument('lyra') | a.Instrument('rhessi'))  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    3 Results from the LYRAClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str6     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2012-03-04 00:00:00 2012-03-06 00:00:00 Proba2       lyra        nan
    2012-03-04 00:00:00 2012-03-06 00:00:00 Proba2       lyra        nan
    2012-03-04 00:00:00 2012-03-06 00:00:00 Proba2       lyra        nan
    <BLANKLINE>
    3 Results from the RHESSIClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str6     str6       str3
    ------------------- ------------------- ------ ---------- ----------
    2012-03-04 00:00:00 2012-03-04 23:59:59 rhessi     rhessi        nan
    2012-03-05 00:00:00 2012-03-05 23:59:59 rhessi     rhessi        nan
    2012-03-06 00:00:00 2012-03-06 23:59:59 rhessi     rhessi        nan
    <BLANKLINE>
    <BLANKLINE>


Indexing search results
-----------------------

The `~sunpy.net.fido_factory.UnifiedResponse` that Fido returns can be
indexed to access a subset of the search results. When doing this, the
results should be treated as a two-dimensional array in which the first
dimension corresponds to the clients which have returned results and the
second to the records returned.

For example, the following code returns a response containing LYRA data from
the `~sunpy.net.dataretriever.sources.LYRAClient`, and EVE data from the
`~sunpy.net.vso.VSOClient`::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"),
    ...                       a.Instrument("lyra") | a.Instrument("eve"))  # doctest: +REMOTE_DATA

If you then wanted to inspect just the LYRA data for the whole time range
specified in the search, you would index this response to see just the
results returned by the `~sunpy.net.dataretriever.sources.LYRAClient`::

    >>> results[0, :]  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the LYRAClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str6     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2012-01-01 00:00:00 2012-01-02 00:00:00 Proba2       lyra        nan
    2012-01-01 00:00:00 2012-01-02 00:00:00 Proba2       lyra        nan
    <BLANKLINE>
    <BLANKLINE>

Or, equivalently::

    >>> results[0]  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the LYRAClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str6     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2012-01-01 00:00:00 2012-01-02 00:00:00 Proba2       lyra        nan
    2012-01-01 00:00:00 2012-01-02 00:00:00 Proba2       lyra        nan
    <BLANKLINE>
    <BLANKLINE>


Normal slicing operations work as with any other Python sequence, e.g.
``results[1,::10]`` to access every tenth file in the result returned by
the second client.

Note that the first (client) index is still necessary even if results
are only found for a single client. So in this case the first result
would be ``results[0,0]`` rather than ``results[0]`` (the latter would return
all results from the first - and only - client and is therefore the
same as ``results``).

.. _downloading_data:

Downloading data
----------------
Once you have located your files via a
`Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>`, you can
download them via `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`::

    >>> downloaded_files = Fido.fetch(results)  # doctest: +SKIP

This downloads the files to the location set in you sunpy config
file.  It also returns a list ``downloaded_files``, of absolute file paths
of where the files have been downloaded to.

You can also specify the path to which you want the data downloaded::

  >>> downloaded_files = Fido.fetch(results, path='/ThisIs/MyPath/to/Data/{file}.fits')  # doctest: +SKIP

This downloads the query results into the directory
``/ThisIs/MyPath/to/Data``, naming each downloaded file with the
filename ``{file}`` obtained from the client, and appended with the suffix
``.fits``. You can also use other properties of the returned query
to define the path where the data is saved.  For example, to save the
data to a subdirectory named after the instrument, use

    >>> downloaded_files = Fido.fetch(results, path='./{instrument}/{file}.fits')  # doctest: +SKIP

You can see the list of options that can be specified in path for all the files
to be downloaded with ``results.response_block_properties``.

.. _VSO: https://sdac.virtualsolar.org/cgi/search
.. _JSOC: http://jsoc.stanford.edu/

------------------------
Using SunPy's HEK module
------------------------

The Heliophysics Event Knowledgebase (HEK) is a repository of feature
and event information about the Sun.  Entries are generated both
by automated algorithms and human observers.  SunPy accesses this
information through the `hek` module, which was developed through
support from the European Space Agency Summer of Code in Space
(ESA-SOCIS) 2011.

1. Setting up the client
------------------------

SunPy's HEK module is in sunpy.net.  It can be imported into your
session as follows:

    >>> from sunpy.net import hek
    >>> client = hek.HEKClient()  # doctest: +REMOTE_DATA

This creates a client that we will use to interact with the HEK.

2. A simple query
-----------------

To search the HEK, you need a start time, an end time, and an event
type.  Times are specified as strings or Python datetime objects.  Event types
are specified as upper case, two letter strings, and are identical to
the two letter abbreviations found at the HEK website,
http://www.lmsal.com/hek/VOEvent_Spec.html.

    >>> tstart = '2011/08/09 07:23:56'
    >>> tend = '2011/08/09 12:40:29'
    >>> event_type = 'FL'
    >>> result = client.search(hek.attrs.Time(tstart,tend),hek.attrs.EventType(event_type))  # doctest: +REMOTE_DATA

The first line defines the search start and end times.  The
second line specifies the event type, in this 'FL' or flare.  Line 4
goes out to the web, contacts the HEK, and queries it for the
information you have requested.  Event data for ALL flares available
in the HEK within the time range 2011/08/09 07:23: 56 UT - 2011/08/09
12:40:20 UT will be returned, regardless of which feature recognition
method used to detect the flare.

Let's break down the arguments of client.search.  The first argument::

    hek.attrs.Time(tstart,tend)

sets the start and end times for the query.  The second argument::

    hek.attrs.EventType(event_type)

sets the type of event to look for.  Since we have defined event_type
= 'FL', this sets the query to look for flares.  We could have also
set the flare event type using the syntax::

    hek.attrs.FL

There is more on the attributes of hek.attrs in section 4 of this
guide.


3. The result
-------------

So, how many flare detections did the query turn up?

    >>> len(result)  # doctest: +REMOTE_DATA
    19

The object returned by the above query is a list of Python dictionary
objects. Each dictionary consists of key-value pairs that exactly
correspond to the parameters listed at
http://www.lmsal.com/hek/VOEvent_Spec.html. You can inspect all the
dictionary keys very simply:

    >>> result[0].keys() # doctest:+SKIP
    [u'skel_startc1',
     u'concept',
     u'frm_versionnumber',
     u'hrc_coord',
     u'refs_orig',....

and so on. Remember, the HEK query we made returns all the flares in
the time-range stored in the HEK, regardless of the feature
recognition method.  The HEK parameter which stores the the feature
recognition method is called "frm_name". Using list comprehensions
(which are very cool), it is easy to get a list of the feature
recognition methods used to find each of the flares in the result
object, for example:

    >>> [elem["frm_name"] for elem in result]  # doctest:+SKIP
    [u'asainz',
     u'asainz',
     u'asainz',
     u'asainz',
     u'asainz',
     u'asainz',
     u'asainz',
     u'SSW Latest Events',
     u'SEC standard',
     u'Flare Detective - Trigger Module',
     u'Flare Detective - Trigger Module',
     u'SSW Latest Events',
     u'SEC standard',
     u'Flare Detective - Trigger Module',
     u'Flare Detective - Trigger Module',
     u'Flare Detective - Trigger Module',
     u'Flare Detective - Trigger Module',
     u'Flare Detective - Trigger Module']

It is likely each flare on the Sun was actually detected multiple
times by many different methods.

4. More complex queries
-----------------------

The HEK client allows you to make more complex queries.  There are two
key features you need to know in order to make use of the full power
of the HEK client.  Firstly, the attribute module - hek.attrs -
describes ALL the parameters stored by the HEK as listed in
http://www.lmsal.com/hek/VOEvent_Spec.html, and the HEK client makes
these parameters searchable.

To explain this, let's have a closer look at hek.attrs. The help
command is your friend here; scroll down to section DATA you will see:

    >>> help(hek.attrs) # doctest:+SKIP
    AR = <sunpy.net.hek.attrs.AR object>
    Area = <sunpy.net.hek.attrs.Area object>
    Bound = <sunpy.net.hek.attrs.Bound object>
    BoundBox = <sunpy.net.hek.attrs.BoundBox object>
    CC = <sunpy.net.hek.attrs.CC object>
    CD = <sunpy.net.hek.attrs.CD object>
    CE = <sunpy.net.hek.attrs.CE object>
    CH = <sunpy.net.hek.attrs.EventType object>
    CJ = <sunpy.net.hek.attrs.EventType object>
    CR = <sunpy.net.hek.attrs.EventType object>
    CW = <sunpy.net.hek.attrs.EventType object>
    EF = <sunpy.net.hek.attrs.EF object>
    ER = <sunpy.net.hek.attrs.EventType object>
    Event = <sunpy.net.hek.attrs.Event object>
    FA = <sunpy.net.hek.attrs.EventType object>
    FE = <sunpy.net.hek.attrs.EventType object>
    FI = <sunpy.net.hek.attrs.FI object>
    FL = <sunpy.net.hek.attrs.FL object>
    FRM = <sunpy.net.hek.attrs.FRM object>
    etc etc...

The object hek.attrs knows the attributes of the HEK.  You'll see that
one of the attributes is a flare object::

    FL = <sunpy.net.hek.attrs.FL object>

We can replace hek.attrs.EventType('FL') with hek.attrs.FL - they do
the same thing, setting the query to look for flare events.  Both
methods of setting the event type are provided as a convenience

Let's look further at the FRM attribute::

    >>> help(hek.attrs.FRM) # doctest:+SKIP
    Help on FRM in module sunpy.net.hek.attrs object:
    class FRM(__builtin__.object)
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  Contact = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  HumanFlag = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  Identifier = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  Institute = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  Name = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  ParamSet = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  SpecificID = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  URL = <sunpy.net.hek.attrs._StringParamAttrWrapper object>
     |
     |  VersionNumber = <sunpy.net.hek.attrs._StringParamAttrWrapper object>

Let's say I am only interested in those flares identified by the SSW
Latest Events tool.  I can retrieve those entries only from the HEK
with the following command:

    >>> result = client.search( hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type), hek.attrs.FRM.Name == 'SSW Latest Events')  # doctest: +REMOTE_DATA
    >>> len(result)  # doctest: +REMOTE_DATA
    2

We can also retrieve all the entries in the time range which were not
made by SSW Latest Events with the following command:

    >>> result = client.search( hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type),hek.attrs.FRM.Name != 'SSW Latest Events')  # doctest: +REMOTE_DATA
    >>> len(result)  # doctest: +REMOTE_DATA
    17

We are using Python's comparison operators to filter the returns from
the HEK client.  Other comparisons are possible.  For example, let's
say I want all the flares that have a peak flux of over 4000.0:

    >>> result = client.search(hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type), hek.attrs.FL.PeakFlux > 4000.0)  # doctest: +REMOTE_DATA
    >>> len(result)  # doctest: +REMOTE_DATA
    1

Multiple comparisons can be included.  For example, let's say I want
all the flares with a peak flux above 1000 AND west of 800 arcseconds
from disk center of the Sun

    >>> result = client.search(hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type), hek.attrs.Event.Coord1 > 800, hek.attrs.FL.PeakFlux > 1000.0)  # doctest: +REMOTE_DATA

Multiple comparison operators can be used to filter the results back
from the HEK.

The second important feature about the HEK client is that the
comparisons we've made above can be combined using Python's logical
operators.  This makes complex queries easy to create.  However, some
caution is advisable.  Let's say I want all the flares west of 50
arcseconds OR have a peak flux over 1000.0:

    >>> result = client.search(hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type), (hek.attrs.Event.Coord1 > 50) or (hek.attrs.FL.PeakFlux > 1000.0) )  # doctest: +REMOTE_DATA

and as a check

    >>> [elem["fl_peakflux"] for elem in result] # doctest: +REMOTE_DATA
    [None, None, None, None, None, None, None, 2326.86, 1698.83, None, None, 2360.49, 3242.64, 1375.93, 6275.98, 923.984, 1019.83]

    >>> [elem["event_coord1"] for elem in result] # doctest: +REMOTE_DATA
    [51.0, 51.0, 51.0, 924.0, 924.0, 924.0, 69.0, 883.2, 883.2, 69.0, 69.0, 883.2, 883.2, 883.2, 883.2, 883.2, 883.2]

Note that some of the fluxes are returned as "None".  This is because
some feature recognition methods for flares do not report the peak
flux.  However, because the location of event_coord1 is greater than
50, the entry from the HEK for that flare detection is returned.

Let's say we want all the flares west of 50 arcseconds AND have a peak
flux over 1000.0:

    >>> result = client.search(hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type), (hek.attrs.Event.Coord1 > 50) and (hek.attrs.FL.PeakFlux > 1000.0) )  # doctest: +REMOTE_DATA

    >>> [elem["fl_peakflux"] for elem in result] # doctest: +REMOTE_DATA
    [2326.86, 1698.83, 2360.49, 3242.64, 1375.93, 6275.98, 1019.83]
    >>> [elem["event_coord1"] for elem in result] # doctest: +REMOTE_DATA
    [883.2, 883.2, 883.2, 883.2, 883.2, 883.2, 883.2]

In this case none of the peak fluxes are returned with the value
`None`.  Since we are using an `and` logical operator we need a result
from the `(hek.attrs.FL.PeakFlux > 1000.0)` filter.  Flares that have
`None` for a peak flux cannot provide this, and so are excluded.  The
`None` type in this context effectively means "Don't know"; in such
cases the client returns only those results from the HEK that
definitely satisfy the criteria passed to it.

5. Getting data for your event
------------------------------

The 'hek2vso' module allows you to take an HEK event and acquire VSO
records specific to that event.

    >>> from sunpy.net import hek2vso
    >>> h2v = hek2vso.H2VClient()  # doctest: +REMOTE_DATA

There are several ways to use this capability.  For example, you can
pass in a list of HEK results and get out the corresponding VSO
records.  Here are the VSO records returned via the tenth result from
the HEK query in Section 2 above:

    >>> result = client.search(hek.attrs.Time(tstart,tend),hek.attrs.EventType(event_type))  # doctest: +REMOTE_DATA
    >>> vso_records = h2v.translate_and_query(result[10])  # doctest: +REMOTE_DATA
    >>> len(vso_records[0])  # doctest: +REMOTE_DATA
    31

Result 10 is an HEK entry generated by the "Flare Detective" automated
flare detection algorithm running on the AIA 193 angstrom waveband.
The VSO records are for full disk AIA 193 images between the start and
end times of this event.  The 'translate_and_query' function uses
exactly that information supplied by the HEK in order to find the
relevant data for that event.  Note that the VSO does not generate
records for all solar data, so it is possible that an HEK entry
corresponds to data that is not accessible via the VSO.

You can also go one step further back, passing in a list of HEK
attribute objects to define your search, the results of which are then
used to generate their corresponding VSO records:

   >>> from sunpy.net import hek
   >>> q = h2v.full_query((hek.attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'), hek.attrs.EventType('FL')))  # doctest: +REMOTE_DATA

The full capabilities of the HEK query module can be used in this
function (see above).

Finally, for greater flexibility, it is possible to pass in a list of
HEK results and create the corresponding VSO query attributes.

    >>> vso_query = hek2vso.translate_results_to_query(result[10:11])  # doctest: +REMOTE_DATA
    >>> vso_query[0]  # doctest: +REMOTE_DATA
    [<Time(<Time object: scale='utc' format='isot' value=2011-08-09T07:22:44.000>, <Time object: scale='utc' format='isot' value=2011-08-09T07:28:56.000>, None)>, <Source('SDO')>, <Instrument('AIA')>, <Wavelength(193.0, 193.0, 'Angstrom')>]

This function allows users finer-grained control of VSO queries
generated from HEK results.

The 'hek2vso' module was developed with support from the 2013 Google
Summer of Code.


---------------------------------------
Querying and Downloading Data from JSOC
---------------------------------------

Joint Science Operations Center (JSOC) contains data products from the Solar Dynamics Observatory,
as well as certain other missions and instruments. These data are available from the JSOC database,
which can be directly accessed by the online `JSOC interface <http://jsoc.stanford.edu/ajax/lookdata.html>`_.

SunPy's JSOC Client provides an easier interface to query for JSOC data and make export requests.
It uses `drms module <http://docs.sunpy.org/projects/drms>`_ as its backend, and exposes a similar API as
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
-----

SunPy's Fido module is in `sunpy.net`.  It can be imported as follows:

    >>> from sunpy.net import Fido, attrs as a

The JSOC client handles the particulars of how the data from
the data provider is downloaded to your computer.

.. warning::

    You must have an email address registered with JSOC before you are allowed to make a request.
    See `this <http://jsoc.stanford.edu/ajax/register_email.html>` to register your email address.

Querying the JSOC
-----------------

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
^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's start with a very simple query.  We could ask for all ``hmi.v_45s`` series data
between January 1st from 00:00 to 01:00, 2014.
We can add email address to the search query with the :mod:`jsoc.Notify` attribute.
Please note you can search without this but right now, you can not add the email address after the search.

    >>> res = Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Series('hmi.v_45s'),
    ...                   a.jsoc.Notify('sunpy@sunpy.org'))  # doctest: +REMOTE_DATA

This returns an `~sunpy.net.fido_factory.UnifiedResponse` object containing
information on the available online files which fit the criteria specified by
the attrs objects in the above call. It does not download the files.

To see a summary of results of our query, simply type the name of the
variable set to the Fido search, in this case, res::

    >>> res  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    81 Results from the JSOCClient:
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
             str23            str7     str10    float64   int64
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

    >>> a.Time(astropy.time.Time('2014-01-01T00:00:00', scale='tai'), astropy.time.Time('2014-01-01T01:00:00', scale='tai'))  # doctest: +SKIP

The second argument::

    >>> a.jsoc.Series('hmi.v_45s')  # doctest: +SKIP

sets the series we are looking for.

So what is going on here?
The notion is that a JSOC query has a set of attribute objects, imported as ``a.jsoc``,
that are specified to construct the query.

``a.jsoc.Series()`` is compulsory to be provided in each of the jsoc queries. Apart from this,
at least one PrimeKey must be passed (generally ``a.Time()``).

The third argument::

    >>> a.jsoc.Notify('sunpy@sunpy.org')  # doctest: +SKIP

tells JSOC what email address you are registered with and to email when your request is ready to download.

Querying with other PrimeKeys
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Other than Time, one other PrimeKey is supported with in-built attribute.
In case of AIA series, ``a.jsoc.Wavelength()`` can be passed as a PrimeKey::

    >>> import astropy.units as u
    >>> res = Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                               a.jsoc.Notify('sunpy@sunpy.org'),
    ...                               a.jsoc.Series('aia.lev1_euv_12s'),
    ...                               a.jsoc.Wavelength(304*u.AA))  # doctest: +REMOTE_DATA

Note that, only Time and Wavelength are in-built attributes here. If you need to pass any other PrimeKey,
it should be passed like this::

    >>> a.jsoc.PrimeKey('HARPNUM', '4864')  # doctest: +SKIP

If 2 or more PrimeKeys need to be passed together::

    >>> a.jsoc.PrimeKey('HARPNUM', '4864') & a.jsoc.PrimeKey('CAMERA', '2')  # doctest: +SKIP

Also, note that the pre-defined primekeys, Time and Wavelength can also be passed as above, but you need to
specify the exact keyword for it. For e.g. by::

    >>> a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.PrimeKey('WAVELNTH', '161')  # doctest: +SKIP

If the correct keyword is not specified, or the passed PrimeKey is not supported by the given series, a
meaningful error will be thrown, which will give you the PrimeKeys supported by that series. Hence, by looking
at the error, one can easily retry building the query with correct PrimeKeys.

Another important thing to note is that, Wavelength when passed through in-built attribute, should be passed as an
Astropy quantity. Specifying spectral units in arguments is necessary or an error will be raised.
For more information on units, see `~astropy.units`.
But, when the same is passed through PrimeKey attribute, it should be passed as a string. All
other PrimeKey values passed through PrimeKey attribute, must be passed as a string.


Manually specifying keyword data to fetch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Upon doing ``Fido.search()`` as described above, only a limited set of keywords are returned in the response
object. These default keywords are ``'DATE'``, ``'TELESCOP'``, ``'INSTRUME'``, ``'T_OBS'`` and ``'WAVELNTH'``.

If you want to get a manual set of keywords in the response object, you can pass the set of keywords using
`~sunpy.net.jsoc.attrs.Keys` attribute.

    >>> res = Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                   a.jsoc.Series('hmi.v_45s'), a.jsoc.Notify('sunpy@sunpy.org'),
    ...                   a.jsoc.Keys(['TELESCOP', 'INSTRUME', 'T_OBS']))  # doctest: +REMOTE_DATA

The parameter passed into ``a.jsoc.Keys()`` can be either a list of strings, or a string with keywords seperated by
comma and a space. Meaning to say,: ``a.jsoc.Keys(['TELESCOP', 'INSTRUME', 'T_OBS'])`` and
``jsoc.attrs.Keys('TELESCOP, INSTRUME, T_OBS')``

both are correct.

Passing an incorrect keyword won't throw an error, but the corresponding column in the table will
contain ``Invalid KeyLink``.

To get all of the keywords, you can either use the `~sunpy.net.jsoc.JSOCClient.search_metadata` method,
or alternatively pass ``a.jsoc.Keys('***ALL***')`` with the series name and PrimeKey.


Using Segments
^^^^^^^^^^^^^^
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
    ...             a.jsoc.Segment('continuum') & a.jsoc.Segment('magnetogram'))  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    61 Results from the JSOCClient:
             T_REC          TELESCOP  INSTRUME WAVELNTH CAR_ROT
             str23            str7      str9   float64   int64
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
    <BLANKLINE>
    <BLANKLINE>


Using Sample
^^^^^^^^^^^^
In case you need to query for data, at some interval of time, say every 10 min, you can pass it
using `~sunpy.net.attrs.Sample`. In other words, if you need to query for ``hmi.v_45s`` series data
between January 1st from 00:00 to 01:00, 2014, every 10 minutes, you can do::

    >>> import astropy.units as u
    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Notify('sunpy@sunpy.org'),
    ...             a.jsoc.Series('hmi.v_45s'), a.Sample(10*u.min))  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    7 Results from the JSOCClient:
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
             str23            str7     str10    float64   int64
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Complex queries can be built using ``OR`` operators.

Let's look for 2 different series data at the same time::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'), a.jsoc.Notify('sunpy@sunpy.org'),
    ...             a.jsoc.Series('hmi.v_45s') | a.jsoc.Series('aia.lev1_euv_12s'))  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    81 Results from the JSOCClient:
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
             str23            str7     str10    float64   int64
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
    <BLANKLINE>
    2107 Results from the JSOCClient:
           T_REC         TELESCOP INSTRUME WAVELNTH CAR_ROT
           str20           str7     str5    int64    int64
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
    <BLANKLINE>
    <BLANKLINE>

The two series names are joined together by the operator ``|``.
This is the ``OR`` operator.  Think of the above query as setting a set
of conditions which get passed to the JSOC.  Let's say you want all the
``hmi.v_45s`` data from two separate days::

    >>> Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00') |
    ...             a.Time('2014-01-02T00:00:00', '2014-01-02T01:00:00'),
    ...             a.jsoc.Series('hmi.v_45s'), a.jsoc.Notify('sunpy@sunpy.org'))  # doctest: +REMOTE_DATA +ELLIPSIS
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    81 Results from the JSOCClient:
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
             str23            str7     str10    float64   int64
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
    <BLANKLINE>
    81 Results from the JSOCClient:
             T_REC          TELESCOP  INSTRUME  WAVELNTH CAR_ROT
             str23            str7     str10    float64   int64
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
    <BLANKLINE>
    <BLANKLINE>

Each of the arguments in this query style can be thought of as
setting conditions that the returned records must satisfy.

It should be noted that ``AND`` operator is supported by some of the attributes only. The attributes which
support "&" are `~sunpy.net.jsoc.attrs.PrimeKey` and `~sunpy.net.jsoc.attrs.Segment`.
Using "&" with any other attributes will throw an error.

Downloading data
----------------

To download the files located by `~sunpy.net.fido_factory.UnifiedDownloaderFactory.search`,
you can download them by `~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`::

    >>> downloaded_files = Fido.fetch(res)  # doctest: +SKIP

Using JSOCClient for complex usage
----------------------------------

Fido interface uses `~sunpy.net.jsoc.JSOCClient` in its backend, and combines
the last 2 stages the JSOC process into one. You can directly use the JSOC
client to make queries, instead of the Fido client. This will allow you to
separate the 3 stages of the JSOC process, and perform it individually, hence
allowing a greater control over the whole process.

Setup
^^^^^

SunPy's JSOC module is in `~sunpy.net`.  It can be imported as follows::

    >>> from sunpy.net import jsoc
    >>> client = jsoc.JSOCClient()  # doctest: +REMOTE_DATA

This creates your client object.


Making a query
^^^^^^^^^^^^^^

Querying JSOC using the JSOC client is very similar to what we were doing with Fido.
As above, we have to make sure we have an email address registered with JSOC before you are allowed to make a request.
See `this <http://jsoc.stanford.edu/ajax/register_email.html>`__ to register your email address.
We can add an email address to the search query with the :mod:`jsoc.Notify` attribute.
Please note you can search without this but right now, you can not add the email address after the search::

    >>> from sunpy.net import attrs as a
    >>> res = client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
    ...                     a.jsoc.Series('hmi.v_45s'),
    ...                     a.jsoc.Notify('sunpy@sunpy.org'))  # doctest: +REMOTE_DATA

Apart from the function name, everything is the same. You need to pass the same values in the
`~sunpy.net.jsoc.JSOCClient.search` as you did in `~sunpy.net.fido_factory.UnifiedDownloaderFactory.search`.
Complex queries can be built in a similar way, and all other things are the same.

Staging the request
^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^

Once the status code is 0 you can download the data using the
`~sunpy.net.jsoc.JSOCClient.get_request` method::

    >>> res = client.get_request(requests)  # doctest: +SKIP

This returns a Results instance which can be used to watch the progress of the download::

    >>> res.wait(progress=True)   # doctest: +SKIP

.. automodapi:: sunpy.net
   :no-heading:

.. automodapi:: sunpy.net.fido_factory


Dataretriever
-------------

.. automodapi:: sunpy.net.dataretriever
   :allowed-package-names: sources
   :headings: ^#

.. automodapi:: sunpy.net.dataretriever.sources
   :headings: #~


VSO
---

.. automodapi:: sunpy.net.vso
   :headings: ^#

.. automodapi:: sunpy.net.vso.attrs
   :headings: #~

JSOC
----

.. automodapi:: sunpy.net.jsoc
    :headings: ^#

.. automodapi:: sunpy.net.jsoc.attrs
    :headings: #~


HEK
---

.. automodapi:: sunpy.net.hek
    :headings: ^#

.. automodapi:: sunpy.net.hek2vso
    :headings: ^#


HELIO
-----

.. automodapi:: sunpy.net.helio
    :headings: ^#

.. automodapi:: sunpy.net.helio.hec
    :headings: #~

Helioviewer
-----------

.. automodapi:: sunpy.net.helioviewer
    :headings: ^#

Attr
----

.. automodapi:: sunpy.net.attr
   :no-heading:
