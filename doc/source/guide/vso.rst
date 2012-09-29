------------------------
Using SunPy's VSO module
------------------------

The Virtual Solar Observatpry (VSO) is a service which presents a
homogenoeous interface to heterogeneous data-sets and services.  Using
the VSO, a user can query multiple data providers simultaneously, and
then download the relevant data.  SunPy uses the VSO through the 'vso'
module, which was developed through support from the European Space
Agency Summer of Code in Space (ESA-SOCIS) 2011.

1. Setting up the VSO interaction
--------------------------

SunPy's VSO module is in sunpy.net.  It can be imported into your
IPython session as follows:

    >>> from sunpy.net import vso

We will also be using the datetime object from the datetime module, so
we'll import that too

    >>> from datetime import datetime



2. A simple query - using the legacy query
---------------------------------

Obtaining data via the VSO is essentially a two-stage process.  In the
first stage, you ask the VSO to find the data you want.  The VSO
queries various data-providers looking for your data.  In the second
stage, you download the data, if there is any data that matches your
request.  The VSO client handles the particulars of how the data from
the data provider is downloaded to your computer.

Let's start with a very simple query.  To search the VSO, you need a
start time, an end time, and an instrument. We have provided two
different syntaxes for doing this search.  The first query syntax
basically copies what you already may be used to from Solarsoft/IDL's
VSO query client, VSO_SEARCH.pro.  This is known as a 'legacy' query,
purely because the syntax is based on the legacy of Solarsoft/IDL's
VSO query client.  The second query syntax is much more powerful, but
slightly less familiar to Solarsoft/IDL users (which is why we have
two different syntaxes).

Let's say I want all the EIT data between 2001/01/01 and 2001/01/02.
Using the legacy query syntax, this is simply

    >>> client = vso.VSOClient()
    >>> qr = client.query_legacy(tstart = '2001/01/01', tend =
    '2001/01/02', instrument = 'EIT')

which is almost identical to what you would type in a Solarsoft/IDL
session.  So, what's happening with this command?  The client is going
out to the web to query the VSO to ask how many files EIT images are
in the archive between the start of 2001/01/01 and the start of
2001/01/02.

The same query can also be performed using a slightly different
syntax.  For example

    >>> qr = client.query_legacy(tstart = datetime(2001,1,1), tend =
    datetime(2001,1,2), instrument = 'EIT')

and 

    >>> qr = client.query_legacy(datetime(2001,1,1),
        datetime(2001,1,2), instrument = 'EIT')

both give the same result.

The variable "qr" is a Python list of response objects, each one
of which is a record found by the VSO. How many records have been
found?  You can find that out be typing

    >>> qr.num_records()
    122

To get a little bit more information, try

    >>> qr.show()


The Solarsoft legacy query has more keywords available: to find out
more about the legacy query, type: 

    >>> help(client.query_legacy)

As an example, let's say you just want the EIT 171 Angstrom files for
that data.  These files can be found by

    >>> qr = client.query_legacy(tstart = '2001/01/01', tend =
    '2001/01/02', instrument = 'EIT', min_wave = '171', max_wave =
    '171', unit_wave = 'Angstrom')

which yields four results, the same as the VSO IDL client.

Having located the data you want, you can download it using the
following command:

    >>>> res = client.get(qr, path = '/Users/ireland/Desktop/Data/{file}.fits')

This downloads the query results into the directory
/Users/ireland/Dekstop/Data naming downloaded file with the filename
'{file}' obtained from the VSO , and appended with the suffix '.fits'.

The '{file}' option uses the file name obtained by the VSO for each
file.  You can also use other properties of the query return to define
the path where the data is saved.  For example,





Using the
legacy query keywords it is very easy to translate a Solarsoft/IDL VSO
command into the equivalent SunPy VSO legacy query.  However, more
powerful queries are possible with the new query style, which is
descibed below.


2. The new query style
------------------



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

    >>> len(result)
    19

The object returned by the above query is a list of Python dictionary
objects.  Each dictionary consists of key-value pairs that exactly
correspond to the parameters listed at
http://www.lmsal.com/hek/VOEvent_Spec.html. You can inspect all the
dictionary keys very simply:

    >>> result[0].keys()
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

    >>> [elem["frm_name"] for elem in result] 
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

    >>> help(hek.attrs)
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

Let's look further at the FRM attribute:

    >>> help(hek.attrs.FRM)
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

    >>> [elem["fl_peakflux"] for elem in result]
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

    >>> [elem["event_coord1"] for elem in result]
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

    >>> [elem["fl_peakflux"] for elem in result] 
    [2326.86, 1698.83, 2360.49, 3242.64, 1375.93, 6275.98]
    >>> [elem["event_coord1"] for elem in result]
    [883.2, 883.2, 883.2, 883.2, 883.2, 883.2]

In this case none of the peak fluxes are returned with the value
"None".  Since we are using an "and" logical operator we need a result
from the "(hek.attrs.FL.PeakFlux > 1000.0)" filter.  Flares that have
"None" for a peak flux cannot provide this, and so are excluded.  The
"None" type in this context effectively means "Don't know"; in such
cases the client returns only those results from the HEK that
definitely satisfy the criteria passed to it. 


