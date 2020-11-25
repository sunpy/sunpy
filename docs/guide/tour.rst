A brief tour of SunPy
*********************

This brief tutorial will walk you through some
of the functionality offered by SunPy. Start by reading this tutorial
and trying out some of the examples demonstrated. Once you've completed the
tutorial check out the rest of the :doc:`User Guide </guide/index>` for a more
thorough look at the functionality available.

Sample Data
===========
This tour makes use of a number of sample data files which you will need to
download. This will happen when the sample data is imported for the first time.

Maps
====
Maps are the primary data type in SunPy. They are spatially aware data arrays.
There are maps for a 2D image, a time series of 2D images or temporally aligned
2D images.

**Creating a Map**

SunPy supports many different data products from various sources 'out of the
box'. We shall use SDO's AIA instrument as an example in this tutorial. The
general way to create a Map from one of the supported data products is with the
`~sunpy.map.Map` function from the `sunpy.map` submodule.
`~sunpy.map.Map` takes either a filename, a list of
filenames or a data array and header. We can test
`~sunpy.map.Map` with:


.. plot::
    :include-source:

    import sunpy.data.sample
    import sunpy.map

    aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    aia.peek()

This returns a map named ``aia`` which can be manipulated with standard SunPy map commands.
For more information about maps checkout the :doc:`map guide <data_types/maps>`
and the :ref:`map`.

TimeSeries
==========

SunPy handles time series data, fundamental to the study of any real world
phenomenon, by creating a TimeSeries object. A timeseries consists of two parts;
times and measurements taken at those times. The data can either be in your
current Python session, alternatively within a local or remote file.
In the code block that follows, data is taken from a file containing samples
from a file containing samples from the GOES satellite's X-ray Sensors (XRS).


.. plot::
    :include-source:

    import numpy as np
    import sunpy.data.sample
    import sunpy.timeseries as ts

    my_timeseries = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')
    my_timeseries.peek()

We've created this timeseries object by passing TimeSeries a string which
represents the name of a GOES lightcurve file. The
`.peek() <sunpy.timeseries.GenericTimeSeries.peek>` method plots the timeseries
data and displays the plot with some default settings. You can also use
`my_timeseries.plot() <sunpy.timeseries.GenericTimeSeries.plot>` if you want more
control over the style of the output plot.

For more information about TimeSeries, check out the
:doc:`timeseries guide <data_types/timeseries>` and the
and the :ref:`timeseries_code_ref`.

Plotting
========

SunPy uses a matplotlib-like interface to its plotting so more complex plots can
be built by combining SunPy with matplotlib. If you're not familiar with
plotting in matplotlib, you should `learn the basics <https://matplotlib.org/users/tutorials.html>`__
before continuing with this guide.

Let's begin by creating a simple plot of an AIA image. To make things easy,
SunPy includes several example files which are used throughout the docs. These
files have names like ``sunpy.data.sample.AIA_171_IMAGE`` and ``sunpy.data.sample.RHESSI_IMAGE``.

Try typing the below example into your interactive Python shell.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample

    aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    aia.peek()

If everything has been configured properly you should see an AIA image with
the default AIA 17.1 colormap, a colorbar on the right-hand side and a title and some
labels.

There is lot going on here, but we will walk you through the example. Briefly,
the first line is importing SunPy, and the second importing the sample data
files. On the third line we create a SunPy Map object which is a spatially-aware
image. On the last line we then plot the `~sunpy.map.Map` object, using the built in 'quick plot'
function `~sunpy.map.GenericMap.peek`.

SunPy uses a matplotlib-like interface to it's plotting so more complex
plots can be built by combining SunPy with matplotlib.

.. plot::
    :include-source:

    import sunpy.map
    import matplotlib.pyplot as plt
    import sunpy.data.sample

    aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    fig = plt.figure()
    ax = plt.subplot(111, projection=aia)

    aia.plot()
    aia.draw_limb()
    aia.draw_grid()
    plt.colorbar()

    plt.show()

For more information check out :ref:`plotting`.

Solar Physical Constants
========================

SunPy contains a convenient list of solar-related physical constants. Here is
a short bit of code to get you started: ::

    >>> from sunpy.sun import constants as con

    # one astronomical unit (the average distance between the Sun and Earth)
    >>> print(con.au)
      Name   = Astronomical Unit
      Value  = 149597870700.0
      Uncertainty  = 0.0
      Unit  = m
      Reference = IAU 2012 Resolution B2

    # the solar radius
    >>> print(con.radius)
      Name   = Nominal solar radius
      Value  = 695700000.0
      Uncertainty  = 0.0
      Unit  = m
      Reference = IAU 2015 Resolution B 3

Not all constants have a shortcut assigned to them (as above). The rest of the constants
are stored in a dictionary. The following code grabs the dictionary and gets all of the
keys.::

    >>> solar_constants = con.constants
    >>> solar_constants.keys()
    dict_keys(['mass', 'radius', 'luminosity', 'mean distance',
               'perihelion distance', 'aphelion distance', 'age',
               'solar flux unit', 'visual magnitude', 'average angular size',
               'surface area', 'average density', 'surface gravity',
               'moment of inertia', 'volume', 'escape velocity', 'oblateness',
               'metallicity', 'sunspot cycle', 'average intensity',
               'effective temperature', 'mass conversion rate', 'center density',
               'center temperature', 'absolute magnitude', 'mean energy production',
               'ellipticity', 'GM', 'W_0', 'sidereal rotation rate',
               'first Carrington rotation (JD TT)',
               'mean synodic period', 'alpha_0',
               'delta_0'])

You can also use the function `sunpy.sun.constants.print_all()` to print out a table of all of the values
available. These constants are provided as a convenience so that everyone is using the same
(accepted) values. For more information check out :ref:`sun_code_ref`.

Quantities and Units
====================

Many capabilities in SunPy make use of physical quantities that are specified
with units. SunPy uses `~astropy.units` to implement this functionality.
Quantities and units are powerful tools for keeping track of variables with
physical meaning and make it straightforward to convert the same physical
quantity into different units. To learn more about the capabilities of
quantities and units, consult :ref:`units-coordinates-sunpy` or
`the astropy tutorial <http://learn.astropy.org/Quantities.html>`__.

To demonstrate this, let's look at the solar radius constant. This is a physical quantity
that can be expressed in length units ::

    >>> from sunpy.sun import constants as con
    >>> con.radius
    <<class 'astropy.constants.iau2015.IAU2015'> name='Nominal solar radius' value=695700000.0 uncertainty=0.0 unit='m' reference='IAU 2015 Resolution B 3'>

shows the solar radius in units of meters.  The same physical quantity can be expressed in different units instead using the ``.to()`` method::

    >>> con.radius.to('km')
    <Quantity 695700. km>

or equivalently::

    >>> import astropy.units as u
    >>> con.radius.to(u.km)
    <Quantity 695700. km>

If, as is sometimes the case, you need just the raw value or the unit from a quantity, you can access these individually
with the ``value`` and ```unit`` attributes, respectively::

    >>> r = con.radius.to(u.km)
    >>> r.value
    695700.0
    >>> r.unit
    Unit("km")

This is useful, but the real power of units is in using them in calculations.
Suppose you have the radius of a circle and would like to calculate its area.
The following code implements this::

    >>> import numpy as np
    >>> import astropy.units as u

    >>> def circle_area(radius):
    ...     return np.pi * radius ** 2

The first line imports numpy, and the second line imports astropy's units
module. The function then calculates the area based on a given radius. When
it does this, it tracks the units of the input and propagates them through
the calculation. Therefore, if we define the radius in meters, the area will
be in meters squared::

    >>> circle_area(4 * u.m)
    <Quantity 50.26548246 m2>

This also works with different units, for example ::

    >>> circle_area(4 * u.imperial.foot)
    <Quantity 50.26548246 ft2>

As demonstrated above, we can convert between different systems of measurement.
For example, if you want the area of a circle in square feet, but were given
the radius in meters, then you can convert it before passing it into the function::

    >>> circle_area((4 * u.m).to(u.imperial.foot))
    <Quantity 541.05315022 ft2>

or you can convert the output::

    >>> circle_area(4 * u.m).to(u.imperial.foot ** 2)
    <Quantity 541.05315022 ft2>


This is an extremely brief summary of the powerful capbilities of Astropy units.  To find out more, see
the `the astropy tutorial <http://learn.astropy.org/Quantities.html>`__ and
`documentation <https://docs.astropy.org/en/stable/units/index.html>`__


Working with Times
==================

SunPy also contains a number of convenience functions for working with dates
and times. Here is a short example: ::

    >>> import sunpy.time

    # parsing a standard time strings
    >>> sunpy.time.parse_time('2004/02/05 12:00')
    <Time object: scale='utc' format='isot' value=2004-02-05T12:00:00.000>

    # This returns a astropy.time.Time object. All SunPy functions which require
    # time as an input sanitize the input using parse_time.

    # the julian day
    >>> sunpy.time.parse_time((2010,4,30)).jd
    2455316.5

    # TimeRange objects are useful for representing ranges of time
    >>> time_range = sunpy.time.TimeRange('2010/03/04 00:10', '2010/03/04 00:20')
    >>> time_range.center
    <Time object: scale='utc' format='isot' value=2010-03-04T00:15:00.000>

For more information about working with time in SunPy checkout the :doc:`time guide <time>`.


Obtaining Data
==============

SunPy supports searching for and fetching data from a variety of sources,
including the `VSO <https://virtualsolar.org/>`__ and the
`JSOC <http://jsoc.stanford.edu/>`__. The majority of SunPy's clients can be
queried using the `sunpy.net.Fido` interface. An example of searching the VSO using this
is below::

  >>> from sunpy.net import Fido, attrs as a

  >>> results = Fido.search(a.Time("2011-09-20T01:00:00", "2011-09-20T02:00:00"),
  ...                       a.Instrument.eit)   # doctest:  +REMOTE_DATA
  >>> Fido.fetch(results, path="./directory/")  # doctest: +SKIP
  ['./directory/efz20110920.010015',
   './directory/efz20110920.010613',
   './directory/efz20110920.011353',
   './directory/efz20110920.011947']

For more information and examples of downloading data with SunPy see :ref:`acquiring_data`.

Database Package
================

The database package can be used to keep a local record of all files downloaded
from the VSO, this means that two searches of the VSO which overlap will not
re-download data.

A simple example of this is shown below::


    >>> import astropy.units as u
    >>> from sunpy.net import Fido, attrs as a
    >>> from sunpy.database import Database

    >>> db = Database()
    >>> db.fetch(a.Time("2011-09-20T01:00:00", "2011-09-20T02:00:00"),
    ...          a.Instrument.aia, a.Sample(45*u.min))  # doctest: +REMOTE_DATA
    >>> db.commit()  # doctest: +REMOTE_DATA
    >>> db  # doctest: +SKIP
    <Table length=4>
     id  observation_time_start observation_time_end ...    download_time      size
    str1         str19                 str19         ...        str19          str7
    ---- ---------------------- -------------------- ... ------------------- -------
       1    2011-09-20 01:00:00  2011-09-20 01:00:01 ... 2020-11-21 14:15:30 66200.0
       2    2011-09-20 01:00:00  2011-09-20 01:00:01 ... 2020-11-21 14:15:30 66200.0
       3    2011-09-20 01:45:00  2011-09-20 01:45:01 ... 2020-11-21 14:15:30 66200.0
       4    2011-09-20 01:45:00  2011-09-20 01:45:01 ... 2020-11-21 14:15:30 66200.0

If you then do a second query::

    >>> db.fetch(a.Time("2011-09-20T01:00:00", "2011-09-20T02:45:00"),
    ...          a.Instrument.aia, a.Sample(45*u.min))  # doctest: +REMOTE_DATA
    >>> db.commit()  # doctest: +REMOTE_DATA
    >>> db  # doctest: +SKIP
    <Table length=6>
     id  observation_time_start observation_time_end ...    download_time      size
    str1         str19                 str19         ...        str19          str7
    ---- ---------------------- -------------------- ... ------------------- -------
       1    2011-09-20 01:00:00  2011-09-20 01:00:01 ... 2020-11-21 14:15:30 66200.0
       2    2011-09-20 01:00:00  2011-09-20 01:00:01 ... 2020-11-21 14:15:30 66200.0
       3    2011-09-20 01:45:00  2011-09-20 01:45:01 ... 2020-11-21 14:15:30 66200.0
       4    2011-09-20 01:45:00  2011-09-20 01:45:01 ... 2020-11-21 14:15:30 66200.0
       5    2011-09-20 02:30:00  2011-09-20 02:30:01 ... 2020-11-21 14:17:51 66200.0
       6    2011-09-20 02:30:00  2011-09-20 02:30:01 ... 2020-11-21 14:17:51 66200.

A query can then be performed against the database to get the records::

    >>> entries = db.search(a.Time("2011-09-20T01:45:00", "2011-09-20T02:15:00"), a.Instrument.aia)  # doctest: +REMOTE_DATA
    >>> len(entries)  # doctest: +SKIP
    4

You can see that only two extra records were added to the database.
For more information check out the :ref:`database_guide`.
