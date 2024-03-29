.. _sunpy-tutorial-acquiring-data-hek:

**********************************************
Searching the Heliophysics Event Knowledgebase
**********************************************

The Heliophysics Event Knowledgebase (HEK) is a repository of feature and event information about the Sun.
Entries are generated both by automated algorithms and human observers.
sunpy accesses this information through the `~sunpy.net.hek` module, which was developed through support from the European Space Agency Summer of Code in Space (ESA-SOCIS) 2011.

A simple query
**************

To search the HEK, you need a start time, an end time, and an event type.
Times are specified as strings or Python datetime objects.
Event types are specified as upper case, two letter strings, and are identical to the two letter abbreviations found at the HEK website, http://www.lmsal.com/hek/VOEvent_Spec.html.

.. code-block:: python

   >>> from sunpy.net import attrs as a
   >>> from sunpy.net import Fido

   >>> tstart = '2011/08/09 07:23:56'
   >>> tend = '2011/08/09 12:40:29'
   >>> event_type = 'FL'
   >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type))  # doctest: +REMOTE_DATA

``tstart`` and ``tend`` defines the start and end times of the query, and ``event_type`` specifies the event type which in this example we are searching for flares defined as ``FL``.
``Fido.search`` goes out to the web, contacts the HEK, and queries it for the information you have requested.
Event data for ALL flares available in the HEK within the time range 2011/08/09 07:23:56 UT - 2011/08/09 12:40:20 UT will be returned, regardless of which feature recognition method used to detect the flare.

Let's break down the arguments of ``Fido.search``.
The first argument:

.. code-block:: python

   >>> a.Time(tstart,tend)  # doctest: +SKIP

sets the start and end times for the query.

The second argument:

.. code-block:: python

   >>> a.hek.EventType(event_type)  # doctest: +SKIP

sets the type of event to look for.
Since we have defined ``event_type = 'FL'``, this sets the query to look for flares.
We could have also set the flare event type using the syntax:

.. code-block:: python

   >>> a.hek.FL  # doctest: +SKIP

There is more information on the attributes below.

The result
**********

So, how many flare detections did the query turn up?

The result object returned by ``Fido.search`` is a `~.UnifiedResponse` object which contains all the results from any clients used in the search.
The first thing we need to do is access the results from the HEK client, the only ones for the query we gave:

.. code-block:: python

   >>> len(result['hek'])  # doctest: +REMOTE_DATA
   19

This object is an `astropy.table.Table` object with the columns which correspond to the parameters listed at http://www.lmsal.com/hek/VOEvent_Spec.html.

You can inspect all results very simply:

.. code-block:: python

   >>> result['hek']  # doctest: +SKIP

Remember, the HEK query we made returns all the flares in the time-range stored in the HEK, regardless of the feature recognition method.
The HEK parameter which stores the the feature recognition method is called "frm_name".
We can select just this column:

.. code-block:: python

   >>> result["hek"]["frm_name"]  # doctest: +REMOTE_DATA
   <QueryResponseColumn name='frm_name' dtype='str32' length=19>
                             asainz
                             asainz
                             asainz
                             asainz
                             asainz
                             asainz
                             asainz
                  SSW Latest Events
                               SWPC
   Flare Detective - Trigger Module
   Flare Detective - Trigger Module
                               SWPC
                  SSW Latest Events
   Flare Detective - Trigger Module
   Flare Detective - Trigger Module
   Flare Detective - Trigger Module
   Flare Detective - Trigger Module
   Flare Detective - Trigger Module
   Flare Detective - Trigger Module

It is likely each flare on the Sun was actually detected multiple times by many different methods.

More complex queries
********************

There are two key features you need to know in order to make more complex queries.
Firstly, the attribute module - ``attrs.hek`` - describes all the parameters stored by the HEK as listed in http://www.lmsal.com/hek/VOEvent_Spec.html, and the HEK client makes these parameters searchable.

To explain this, let's have a closer look at ``attrs.hek``.
By using the help command; scroll down to section DATA you will see:

.. code-block:: python

   >>> help(a.hek) # doctest:+REMOTE_DATA
   Help on module sunpy.net.hek.attrs in sunpy.net.hek:
   <BLANKLINE>
   NAME
       sunpy.net.hek.attrs
   <BLANKLINE>
   DESCRIPTION
       Attributes that can be used to construct HEK queries. They are different to
       the VSO ones in that a lot of them are wrappers that conveniently expose
       the comparisons by overloading Python operators. So, e.g., you are able
       to say AR & AR.NumSpots < 5 to find all active regions with less than 5 spots.
       As with the VSO query, you can use the fundamental logic operators AND and OR
       to construct queries of almost arbitrary complexity. Note that complex queries
       result in multiple requests to the server which might make them less efficient.
   <BLANKLINE>
   CLASSES
   ...

You'll see that one of the attributes is a flare object:

.. code-block:: python

    FL = <sunpy.net.hek.attrs.FL object>

We can replace ``a.hek.EventType('FL')`` with ``a.hek.FL`` - they do the same thing, setting the query to look for flare events.
Both methods of setting the event type are provided as a convenience

Let's look further at the FRM attribute:

.. code-block:: python

   >>> help(a.hek.FRM) # doctest:+REMOTE_DATA
   Help on FRM in module sunpy.net.hek.attrs object:
   <BLANKLINE>
   class FRM(builtins.object)
    |  Data descriptors defined here:
    |
    |  __dict__
    |      dictionary for instance variables
    |
    |  __weakref__
    |      list of weak references to the object
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
   <BLANKLINE>

Let's say I am only interested in those flares identified by the SSW Latest Events tool.
One can retrieve those entries only from the HEK with the following command:

.. code-block:: python

   >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), a.hek.FRM.Name == 'SSW Latest Events')  # doctest: +REMOTE_DATA
   >>> len(result[0])  # doctest: +REMOTE_DATA
   2

We can also retrieve all the entries in the time range which were not made by SSW Latest Events with the following command:

.. code-block:: python

   >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), a.hek.FRM.Name != 'SSW Latest Events')  # doctest: +REMOTE_DATA
   >>> len(result[0])  # doctest: +REMOTE_DATA
   19

We are using Python's comparison operators to filter the returns from Fido.
Other comparisons are possible.
For example, let's say I want all the flares that have a peak flux of over 4000.0:

.. code-block:: python

   >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), a.hek.FL.PeakFlux > 4000.0)  # doctest: +REMOTE_DATA
   >>> len(result[0])  # doctest: +REMOTE_DATA
   1

Multiple comparisons can be included.
For example, let's say I want all the flares with a peak flux above 1000 AND west of 800 arcseconds from disk center of the Sun:

.. code-block:: python

   >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), a.hek.Event.Coord1 > 800, a.hek.FL.PeakFlux > 1000)  # doctest: +REMOTE_DATA
   >>> len(result[0])  # doctest: +REMOTE_DATA
   7

Multiple comparison operators can be used to filter the results back from the HEK.

The second important feature is that the comparisons we've made above can be combined using Python's logical operators.
This makes complex queries easy to create.
However, some caution is advisable.
Let's say we want all the flares west of 50 arcseconds OR have a peak flux over 1000:

.. code-block:: python

   >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), (a.hek.Event.Coord1 > 50) or (a.hek.FL.PeakFlux > 1000))  # doctest: +REMOTE_DATA

and as a check:

.. code-block:: python

   >>> result["hek"]["fl_peakflux"] # doctest: +REMOTE_DATA
   <QueryResponseColumn name='fl_peakflux' dtype='object' length=17>
      None
      None
      None
      None
      None
      None
      None
   2326.86
   1698.83
      None
      None
   2360.49
   3242.64
   1375.93
   6275.98
   923.984
   1019.83

   >>> result["hek"]["event_coord1"] # doctest: +REMOTE_DATA
   <QueryResponseColumn name='event_coord1' dtype='float64' length=17>
    51.0
    51.0
    51.0
   924.0
   924.0
   924.0
    69.0
   883.2
   883.2
    69.0
    69.0
   883.2
   883.2
   883.2
   883.2
   883.2
   883.2

Note that some of the fluxes are returned as "None".
This is because some feature recognition methods for flares do not report the peak flux.
However, because the location of "event_coord1" is greater than 50, the entry from the HEK for that flare detection is returned.

Let's say we want all the flares west of 50 arcseconds AND have a peak flux over 1000:

.. code-block:: python

   >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type), (a.hek.Event.Coord1 > 50) and (a.hek.FL.PeakFlux > 1000))  # doctest: +REMOTE_DATA

   >>> result["hek"]["fl_peakflux"] # doctest: +REMOTE_DATA
   <QueryResponseColumn name='fl_peakflux' dtype='float64' length=7>
   2326.86
   1698.83
   2360.49
   3242.64
   1375.93
   6275.98
   1019.83
   >>> result["hek"]["event_coord1"] # doctest: +REMOTE_DATA
   <QueryResponseColumn name='event_coord1' dtype='float64' length=7>
   883.2
   883.2
   883.2
   883.2
   883.2
   883.2
   883.2

In this case none of the peak fluxes are returned with the value `None`.
Since we are using an ``and`` logical operator we need a result from the ``(a.hek.FL.PeakFlux > 1000)`` filter.
Flares that have `None` for a peak flux cannot provide this, and so are excluded.
The `None` type in this context effectively means "Don't know"; in such cases the client returns only those results from the HEK that definitely satisfy the criteria passed to it.

Getting data for your event
***************************

The `sunpy.net.hek2vso` module allows you to take an HEK event and acquire VSO records specific to that event.

.. code-block:: python

   >>> from sunpy.net import hek2vso

   >>> h2v = hek2vso.H2VClient()  # doctest: +REMOTE_DATA

There are several ways to use this capability.
For example, you can pass in a list of HEK results and get out the corresponding VSO records.
Here are the VSO records returned via the tenth result from the HEK query in Section 2 above:

.. code-block:: python

   >>> result = Fido.search(a.Time(tstart,tend), a.hek.EventType(event_type))  # doctest: +REMOTE_DATA
   >>> vso_records = h2v.translate_and_query(result[0][10])  # doctest: +REMOTE_DATA
   >>> len(vso_records[0])  # doctest: +REMOTE_DATA
   31

``result[0][10]`` is the HEK entry generated by the "Flare Detective" automated flare detection algorithm running on the AIA 193 angstrom waveband.
The VSO records are for full disk AIA 193 angstrom images between the start and end times of this event.
The `~sunpy.net.hek2vso.H2VClient.translate_and_query` function uses exactly that information supplied by the HEK in order to find the relevant data for that event.
Note that the VSO does not generate records for all solar data, so it is possible that an HEK entry corresponds to data that is not accessible via the VSO.

You can also go one step further back, passing in a list of HEK attribute objects to define your search, the results of which are then used to generate their corresponding VSO records:

.. code-block:: python

   >>> vso_query = h2v.full_query((a.Time('2011/08/09 07:00:00', '2011/08/09 07:15:00'), a.hek.EventType('FL')))  # doctest: +REMOTE_DATA

The full capabilities of the HEK query module can be used in this function (see above).

Finally, for greater flexibility, it is possible to pass in a list of HEK results and create the corresponding VSO query attributes.

.. code-block:: python

   >>> vso_query = hek2vso.translate_results_to_query(result[0][10])  # doctest: +REMOTE_DATA
   >>> vso_query[0]  # doctest: +REMOTE_DATA
   [<sunpy.net.attrs.Time(2011-08-09 07:22:44.000, 2011-08-09 07:28:56.000)>, <sunpy.net.attrs.Source(SDO: The Solar Dynamics Observatory.) object at ...>, <sunpy.net.attrs.Instrument(AIA: Atmospheric Imaging Assembly) object at ...>, <sunpy.net.attrs.Wavelength(193.0, 193.0, 'Angstrom')>]

This function allows users finer-grained control of VSO queries generated from HEK results.
