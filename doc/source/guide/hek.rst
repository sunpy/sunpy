------------------------
Using SunPy's HEK module
------------------------

The Heliophysics Event Knowledgebase (HEK) is a repository of feature and
event information concerning the Sun.  Entries are generated both by automated
algorithms and human observers.  SunPy accesses this information through the
'hek' module, which was developed through support from the European Space
Agency Summer of Code in Space (ESA-SOCIS) 2011.

1. Setting up the client
------------------------

SunPy's HEK 

    from sunpy.net import hek
    client = hek.HEKClient()

This creates a client that we will use to interact with the HEK.

2. A simple query
-----------------

To search the HEK, you need a start time, an end time, and an event type.  Times
are specified as Python datetime objects.  Event types are specified as upper
case, two letter strings, and are identical to the two letter abbreviations
found at the HEK website, http://www.lmsal.com/hek/VOEvent_Spec.html.

	import datetime
    tstart = datetime.datetime(2011,8,9,7,23,56)
    tend = datetime.datetime(2011,8,9,12,40,29)
    EventType = 'FL'
    result = client.query(hek.attrs.Time(tstart,tend), hek.attrs.EventType(EventType))

The first line imports the datetime module.  The second and third lines define
the search start and end times.  Line 4 specifies the event type, in this 'FL'
or flare.  Line 5 goes out to the web, contacts the HEK, and queries it
for the information you have requested.  Event data for ALL flares available to
the HEK within the time range will be returned, regardless of which
feature recognition method used to detect the flare.

3. The result
-------------

The object returned by the above query is a list of Python dictionary objects.
Each dictionary consists of key-value pairs that exactly correspond to the
parameters listed at http://www.lmsal.com/hek/VOEvent_Spec.html.