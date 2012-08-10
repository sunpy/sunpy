------------------------
Using SunPy's HEK module
------------------------

The Heliophysics Event Knowledgebase (HEK) is a repository of feature
and event information concerning the Sun.  Entries are generated both
by automated algorithms and human observers.  SunPy accesses this
information through the 'hek' module, which was developed through
support from the European Space Agency Summer of Code in Space
(ESA-SOCIS) 2011.

1. Setting up the client
------------------------

SunPy's HEK module is in sunpy.net.  It can be imported into your
IPython session as follows:

    In [1]: from sunpy.net import hek
    In [2]: client = hek.HEKClient()

This creates a client that we will use to interact with the HEK.

2. A simple query
-----------------

To search the HEK, you need a start time, an end time, and an event
type.  Times are specified as Python datetime objects.  Event types
are specified as upper case, two letter strings, and are identical to
the two letter abbreviations found at the HEK website,
http://www.lmsal.com/hek/VOEvent_Spec.html.

    In [3]: import datetime
    In [4]: tstart = datetime.datetime(2011,8,9,7,23,56)
    In [5]: tend = datetime.datetime(2011,8,9,12,40,29)
    In [6]: event_type = 'FL'
    In [7]: result = client.query(hek.attrs.Time(tstart,tend),
                                hek.attrs.EventType(event_type))

The first line in the block of code above ("In [3]:") imports the
datetime module.  The second and third lines define the search start
and end times.  The fourth line specifies the event type, in this 'FL'
or flare.  Line 5 ("In [7]:") goes out to the web, contacts the HEK,
and queries it for the information you have requested.  Event data for
ALL flares available in the HEK within the time range 2011/08/09
07:23: 56 UT - 2011/08/09 12:40:20 UT will be returned, regardless of
which feature recognition method used to detect the flare.

Let's break down the arguments of client.query.  The first argument:

    hek.attrs.Time(tstart,tend)

sets the start and end times for the query.  The second argument:

    hek.attrs.EventType(event_type)

sets the type of event to look for.  Since we have defined event_type
= 'FL', this sets the query to look for flares.  We could have also
set the flare event type using the syntax

    hek.attrs.FL

There is more on the attributes of hek.attrs in section 4 of this
guide.


3. The result
-------------

So, how many flare detections did the query turn up?

    In [8]: len(result)
    Out[8]: 19

The object returned by the above query is a list of Python dictionary
objects.  Each dictionary consists of key-value pairs that exactly
correspond to the parameters listed at
http://www.lmsal.com/hek/VOEvent_Spec.html. You can inspect all the
dictionary keys very simply:

    In [10]: result[0].keys()
    Out[10]:
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

    In [11]: [elem["frm_name"] for elem in result]
    Out[11]: 
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
command is your friend here:

    In [12]: help(hek.attrs)

If you scroll down to section DATA you will see

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
one of the attributes is a flare object

   FL = <sunpy.net.hek.attrs.FL object>

We can replace hek.attrs.EventType('FL') with hek.attrs.FL - they do
the same thing, setting the query to look for flare events.  Both
methods of setting the event type are provided as a convenience

Let's look further at the FRM attribute....

    In [13]: help(hek.attrs.FRM)

which yields

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

    In [14]: result = client.query( hek.attrs.Time(tstart,tend), 
                            hek.attrs.EventType(event_type),
                            hek.attrs.FRM.Name == 'SSW Latest Events')
    In [15]: len(result)
    Out[15]: 2

We can also retrieve all the entries in the time range which were not
made by SSW Latest Events with the following command:

    In [16]: result = client.query( hek.attrs.Time(tstart,tend), 
                            hek.attrs.EventType(event_type),
                            hek.attrs.FRM.Name != 'SSW Latest Events')
    In [17]: len(result)
    Out[17]: 17

We are using Python's comparison operators to filter the returns from
the HEK client.  Other comparisons are possible.  For example, let's
say I want all the flares that have a peak flux of over 4000.0:

    In [18]: result = client.query(hek.attrs.Time(tstart,tend),
                            hek.attrs.EventType(event_type),
                            hek.attrs.FL.PeakFlux > 4000.0)
    In [19]: len(result)
    Out[19]: 1

Multiple comparisons can be included.  For example, let's say I want
all the flares with a peak flux above 1000 AND west of 800 arcseconds
from disk center of the Sun

    In [20]: result = client.query(hek.attrs.Time(tstart,tend),
                            hek.attrs.EventType(event_type),
                            hek.attrs.Event.Coord1 > 800,
                            hek.attrs.FL.PeakFlux > 1000.0)

So, comparison operators can be used to filter the results back from
the HEK.

The second important feature about the HEK client is that the
comparisons we've made above can be combined using Python's logical
operators.  This makes complex queries easy to create.  However, some
caution is advisable.  Let's say I want all the flares west of 50
arcseconds OR have a peak flux over 1000.0:

    In [21]: result = client.query(hek.attrs.Time(tstart,tend),
                            hek.attrs.EventType(event_type),
                            (hek.attrs.Event.Coord1 > 50) or 
                            (hek.attrs.FL.PeakFlux > 1000.0) )
and as a check

    In [22]: [elem["fl_peakflux"] for elem in result]
    Out[22]: [None,
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

    In [23]: [elem["event_coord1"] for elem in result]
    Out[23]: [51,
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

    In [24]: result = client.query(hek.attrs.Time(tstart,tend),
                            hek.attrs.EventType(event_type),
                            (hek.attrs.Event.Coord1 > 50) and 
                            (hek.attrs.FL.PeakFlux > 1000.0) )

    In [25]: [elem["fl_peakflux"] for elem in result] 
    Out[25]: [2326.86, 1698.83, 2360.49, 3242.64, 1375.93, 6275.98]
    In [26]: [elem["event_coord1"] for elem in result]
    Out[26]: [883.2, 883.2, 883.2, 883.2, 883.2, 883.2]

In this case none of the peak fluxes are returned with the value
"None".  Since we are using an "and" logical operator we need a result
from the "(hek.attrs.FL.PeakFlux > 1000.0)" filter.  Flares that have
"None" for a peak flux cannot provide this, and so are excluded.  The
"None" type in this context effectively means "Don't know"; in such
cases the client returns only those results from the HEK that
definitely satisfy the criteria passed to it. 


