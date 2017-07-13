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
    >>> client = hek.HEKClient()

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
    >>> result = client.query(hek.attrs.Time(tstart,tend),hek.attrs.EventType(event_type))

The first line defines the search start and end times.  The
second line specifies the event type, in this 'FL' or flare.  Line 4
goes out to the web, contacts the HEK, and queries it for the
information you have requested.  Event data for ALL flares available
in the HEK within the time range 2011/08/09 07:23: 56 UT - 2011/08/09
12:40:20 UT will be returned, regardless of which feature recognition
method used to detect the flare.

Let's break down the arguments of client.query.  The first argument::

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

    >>> len(result)
    19

The object returned by the above query is a list of Python dictionary
objects.  Each dictionary consists of key-value pairs that exactly
correspond to the parameters listed at
http://www.lmsal.com/hek/VOEvent_Spec.html. You can inspect all the
dictionary keys very simply:

    >>> result[0].keys() # doctest:+SKIP
    [u'skel_startc1',
     u'concept',
     u'frm_versionnumber',
     u'hrc_coord',
     u'refs_orig',....

and so on.  Remember, the HEK query we made returns all the flares in
the time-range stored in the HEK, regardless of the feature
recognition method.  The HEK parameter which stores the the feature
recognition method is called "frm_name". Using list comprehensions
(which are very cool), it is easy to get a list of the feature
recognition methods used to find each of the flares in the result
object, for example:

    >>> [elem["frm_name"] for elem in result] # doctest:+SKIP
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

    >>> result = client.query( hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type), hek.attrs.FRM.Name == 'SSW Latest Events')
    >>> len(result)
    2

We can also retrieve all the entries in the time range which were not
made by SSW Latest Events with the following command:

    >>> result = client.query( hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type),hek.attrs.FRM.Name != 'SSW Latest Events')
    >>> len(result)
    17

We are using Python's comparison operators to filter the returns from
the HEK client.  Other comparisons are possible.  For example, let's
say I want all the flares that have a peak flux of over 4000.0:

    >>> result = client.query(hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type), hek.attrs.FL.PeakFlux > 4000.0)
    >>> len(result)
    1

Multiple comparisons can be included.  For example, let's say I want
all the flares with a peak flux above 1000 AND west of 800 arcseconds
from disk center of the Sun

    >>> result = client.query(hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type), hek.attrs.Event.Coord1 > 800, hek.attrs.FL.PeakFlux > 1000.0)

Multiple comparison operators can be used to filter the results back
from the HEK.

The second important feature about the HEK client is that the
comparisons we've made above can be combined using Python's logical
operators.  This makes complex queries easy to create.  However, some
caution is advisable.  Let's say I want all the flares west of 50
arcseconds OR have a peak flux over 1000.0:

    >>> result = client.query(hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type), (hek.attrs.Event.Coord1 > 50) or (hek.attrs.FL.PeakFlux > 1000.0) )

and as a check

    >>> [elem["fl_peakflux"] for elem in result] # doctest:+SKIP
    [None,
    None,
    None,
    None,
    None,
    None,
    None,
    2326.86,
    1698.83,
    None,
    None,
    2360.49,
    3242.64,
    1375.93,
    6275.98,
    923.984]

    >>> [elem["event_coord1"] for elem in result] # doctest:+SKIP
    [51,
    51,
    51,
    924,
    924,
    924,
    69,
    883.2,
    883.2,
    69,
    69,
    883.2,
    883.2,
    883.2,
    883.2,
    883.2]

Note that some of the fluxes are returned as "None".  This is because
some feature recognition methods for flares do not report the peak
flux.  However, because the location of event_coord1 is greater than
50, the entry from the HEK for that flare detection is returned.

Let's say we want all the flares west of 50 arcseconds AND have a peak
flux over 1000.0:

    >>> result = client.query(hek.attrs.Time(tstart,tend), hek.attrs.EventType(event_type), (hek.attrs.Event.Coord1 > 50) and (hek.attrs.FL.PeakFlux > 1000.0) )

    >>> [elem["fl_peakflux"] for elem in result] # doctest:+SKIP
    [2326.86, 1698.83, 2360.49, 3242.64, 1375.93, 6275.98]
    >>> [elem["event_coord1"] for elem in result] # doctest:+SKIP
    [883.2, 883.2, 883.2, 883.2, 883.2, 883.2]

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
    >>> h2v = hek2vso.H2VClient()

There are several ways to use this capability.  For example, you can
pass in a list of HEK results and get out the corresponding VSO
records.  Here are the VSO records returned via the tenth result from
the HEK query in Section 2 above:

    >>> result = client.query(hek.attrs.Time(tstart,tend),hek.attrs.EventType(event_type))
    >>> vso_records = h2v.translate_and_query(result[10])
    >>> len(vso_records[0])
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
   >>> q = h2v.full_query((hek.attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'), hek.attrs.EventType('FL')))

The full capabilities of the HEK query module can be used in this
function (see above).

Finally, for greater flexibility, it is possible to pass in a list of
HEK results and create the corresponding VSO query attributes.

    >>> vso_query = hek2vso.translate_results_to_query(result[10:11])
    >>> vso_query[0]
    [<Time(datetime.datetime(2011, 8, 9, 7, 22, 44), datetime.datetime(2011, 8, 9, 7, 28, 56), None)>, <Source(u'SDO')>, <Instrument(u'AIA')>, <Wave(193.0, 193.0, 'Angstrom')>]

This function allows users finer-grained control of VSO queries
generated from HEK results.

The 'hek2vso' module was developed with support from the 2013 Google
Summer of Code.
