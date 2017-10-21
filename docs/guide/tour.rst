A brief tour of SunPy
=====================

This brief tutorial will walk you through some
of the functionality offered by SunPy. Start by reading this tutorial
and trying out some of the examples demonstrated. Once you've completed the
tutorial check out the rest of the :doc:`User Guide </guide/index>` for a more
thorough look at the functionality available.

Sample Data
-----------
This tour makes use of a number of sample data files which you will need to
download. This will happen when the sample data is imported for the first time.

Maps
----
Maps are the primary data type in SunPy. They are spatially aware data arrays.
There are maps for a 2D image, a time series of 2D images or temporally aligned
2D images.

**Creating a Map**

SunPy supports many different data products from various sources 'out of the
box'. We shall use SDO's AIA instrument as an example in this tutorial. The
general way to create a Map from one of the supported data products is with the
`Map <sunpy.map.map_factory.MapFactory>` function from the `sunpy.map` submodule.
`Map <sunpy.map.map_factory.MapFactory` takes either a filename, a list of
filenames or a data array and header. We can test
`Map <sunpy.map.map_factory.MapFactory>` with:


.. plot::
    :include-source:

    import sunpy.data.sample
    import sunpy.map

    aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    aia.peek()

This returns a map named `aia` which can be manipulated with standard SunPy map commands.
For more information about maps checkout the :doc:`map guide <data_types/maps>`
and the :ref:`map`.

TimeSeries
----------

SunPy handles time series data, fundamental to the study of any real world
phenomenon, by creating a TimeSeries object. A timeseries consists of two parts;
times and measurements taken at those times. The data can either be in your
current Python session, alternatively within a local or remote file. Let's
create some fake data and pass it into a timeseries object.

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
`my_timeseries.plot() <sunpy.timeseries.TimeSeries.plot>` if you want more
control over the style of the output plot.

For more information about TimeSeries, check out the
:doc:`timeseries guide <data_types/timeseries>` and the
and the :ref:`timeseries_code_ref`.

Spectra
-------

SunPy has spectral support for instruments which have such a capacity. CALLISTO,
an international network of Solar Radio Spectrometers, is a specific example.

.. plot::
    :include-source:

    import sunpy.data.sample
    from sunpy.spectra.sources.callisto import CallistoSpectrogram

    image = CallistoSpectrogram.read(sunpy.data.sample.CALLISTO_SPECTRUM)
    image.peek()

For more information about spectra, check out the :doc:`spectra guide <data_types/spectra>`
and the :ref:`spectra_code_ref`.

Plotting
--------

SunPy uses a matplotlib-like interface to its plotting so more complex plots can
be built by combining SunPy with matplotlib. If you're not familiar with
plotting in matplotlib, you should `learn the basics <http://matplotlib.org/users/tutorials.html>`__
before continuing with this guide.

Let's begin by creating a simple plot of an AIA image. To make things easy,
SunPy includes several example files which are used throughout the docs. These
files have names like `sunpy.data.sample.AIA_171_IMAGE` and `sunpy.data.sample.RHESSI_IMAGE`.

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
image. On the last line we then plot the `Map <sunpy.map.map_factory.MapFactory>` object, using the built in 'quick plot'
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
    aia.draw_limb()
    plt.colorbar()

    plt.show()

For more information check out :ref:`plotting`.

Solar Physical Constants
------------------------

SunPy contains a convenient list of solar-related physical constants. Here is
a short bit of code to get you started: ::

    >>> from sunpy.sun import constants as con

    # one astronomical unit (the average distance between the Sun and Earth)
    >>> print(con.au)
      Name   = Astronomical Unit
      Value  = 1.495978707e+11
      Error  = 0.0
      Units  = m
      Reference = IAU 2012 Resolution B2

    # the solar radius
    >>> print(con.radius)
      Name   = Solar radius
      Value  = 695508000.0
      Error  = 26000.0
      Units  = m
      Reference = Allen's Astrophysical Quantities 4th Ed.

Not all constants have a shortcut assigned to them (as above). The rest of the constants
are stored in a dictionary. The following code grabs the dictionary and gets all of the
keys.::

    >>> solar_constants = con.constants
    >>> solar_constants.keys()   # doctest: +NORMALIZE_WHITESPACE
    ['solar flux unit', 'surface area', 'average density', 'radius', 'surface
    gravity', 'ellipticity', 'visual magnitude', 'center density', 'average
    angular size', 'absolute magnitude', 'sunspot cycle', 'effective
    temperature', 'aphelion distance', 'mean energy production', 'mass
    conversion rate', 'average intensity', 'volume', 'metallicity', 'moment of
    inertia', 'escape velocity', 'perihelion distance', 'GM', 'oblateness',
    'mean distance', 'age', 'mass', 'luminosity', 'center temperature']

You can also use the function `sunpy.constants.print_all()` to print out a table of all of the values
available. These constants are provided as a convenience so that everyone is using the same
(accepted) values. For more information check out :ref:`sun_code_ref`.

Quantities and Units
--------------------

Many capabilities in SunPy make use of physical quantities that are specified
with units. SunPy uses `~astropy.units` to implement this functionality.
Quantities and units are powerful tools for keeping track of variables with
physical meaning and make it straightforward to convert the same physical
quantity into different units. To learn more about the capabilities of
quantities and units, consult :ref:`units-coordinates-sunpy` or
`the astropy tutorial <http://www.astropy.org/astropy-tutorials/Quantities.html>`__.

To demonstrate this, let's look at the solar radius constant. This is a physical quantity
that can be expressed in length units ::

    >>> from sunpy.sun import constants as con
    >>> con.radius
    <Constant name=u'Solar radius' value=695508000.0 error=26000.0 units='m' reference=u"Allen's Astrophysical Quantities 4th Ed.">

shows the solar radius in units of meters.  The same physical quantity can be expressed in different units instead using the `.to()` method::

    >>> con.radius.to('km')
    >>> <Quantity 695508.0 km>

or equivalently::

    >>> import astropy.units as u
    >>> con.radius.to(u.km)
    >>> <Quantity 695508.0 km>

If, as is sometimes the case, you need just the raw value or the unit from a quantity, you can access these individually
with the `value` and `unit` attributes, respectively::

    >>> r = con.radius.to(u.km)
    >>> r.value
    695508.0
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
    <Quantity 50.26548245743669 m2>

This also works with different units, for example ::

    >>> circle_area(4 * u.imperial.foot)
    <Quantity 50.26548245743669 ft2>

As demonstrated above, we can convert between different systems of measurement.
For example, if you want the area of a circle in square feet, but were given
the radius in meters, then you can convert it before passing it into the function::

    >>> circle_area((4 * u.m).to(u.imperial.foot))
    <Quantity 541.0531502245425 ft2>

or you can convert the output::

    >>> circle_area(4 * u.m).to(u.imperial.foot ** 2)
    <Quantity 541.0531502245426 ft2>


This is an extremely brief summary of the powerful capbilities of Astropy units.  To find out more, see
the `the astropy tutorial <http://www.astropy.org/astropy-tutorials/Quantities.html>`__ and
`documentation <http://docs.astropy.org/en/stable/units/index.html>`__


Working with Times
------------------

SunPy also contains a number of convenience functions for working with dates
and times. Here is a short example: ::

    >>> import sunpy.time

    # parsing a standard time strings
    >>> sunpy.time.parse_time('2004/02/05 12:00')
    datetime.datetime(2004, 2, 5, 12, 0)

    # This returns a datetime object. All SunPy functions which require
    # time as an input sanitize the input using parse_time.
    >>> sunpy.time.day_of_year('2004-Jul-05 12:00:02')
    187.50002314814816

    # the julian day
    >>> sunpy.time.julian_day((2010,4,30))
    2455316.5

    # TimeRange objects are useful for representing ranges of time
    >>> time_range = sunpy.time.TimeRange('2010/03/04 00:10', '2010/03/04 00:20')
    >>> time_range.center
    datetime.datetime(2010, 3, 4, 0, 15)

For more information about working with time in SunPy checkout the :doc:`time guide <time>`.


Obtaining Data
--------------

SunPy supports searching for and fetching data from a variety of sources,
including the `VSO <http://virtualsolar.org/>`__ and the
`JSOC <http://jsoc.stanford.edu/>`__. The majority of SunPy's clients can be
queried using the `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` interface. An example of searching the VSO using this
is below::

  >>> from sunpy.net import Fido, attrs as a

  >>> results = Fido.search(a.Time("2011-09-20T01:00:00", "2011-09-20T02:00:00"),
                            a.Instrument('EIT'))   # doctest: +NORMALIZE_WHITESPACE

  <sunpy.net.fido_factory.UnifiedResponse object at 0x7fe70e6c6160>
  Results from 1 Provider:

  4 Results from the VSOClient:
    Start Time [1]       End Time [1]    Source Instrument   Type   Wavelength [2]
                                                                        Angstrom
        str19               str19         str4     str3      str8      float64
  ------------------- ------------------- ------ ---------- -------- --------------
  2011-09-20 01:00:15 2011-09-20 01:00:27   SOHO        EIT FULLDISK 171.0 .. 171.0
  2011-09-20 01:06:13 2011-09-20 01:08:15   SOHO        EIT FULLDISK 284.0 .. 284.0
  2011-09-20 01:13:53 2011-09-20 01:14:05   SOHO        EIT FULLDISK 195.0 .. 195.0
  2011-09-20 01:19:47 2011-09-20 01:20:19   SOHO        EIT FULLDISK 304.0 .. 304.0

  >>> Fido.fetch(results, path="./directory/")
  ['./directory/efz20110920.010015',
   './directory/efz20110920.010613',
   './directory/efz20110920.011353',
   './directory/efz20110920.011947']

For more information and examples of downloading data with SunPy see :ref:`acquiring_data`.

Database Package
----------------

The database package can be used to keep a local record of all files downloaded
from the VSO, this means that two searches of the VSO which overlap will not
re-download data.

A simple example of this is shown below::


  >>> import astropy.units as u
  >>> from sunpy.net import Fido, attrs as a
  >>> from sunpy.database import Database

  >>> db = Database()
  >>> db.fetch(a.Time("2011-09-20T01:00:00", "2011-09-20T02:00:00"),
  ...          a.Instrument('AIA'), a.vso.Sample(15*u.min))
  >>> db.commit()

  >>> db

  <Table length=10>
  id  observation_time_start observation_time_end instrument source provider  physobs  wavemin wavemax                                      path                                              fileid          tags starred    download_time      size
  str2         str19                 str19            str3     str3    str4      str9     str4    str4                                      str77                                              str24           str3   str2         str19          str7
  ---- ---------------------- -------------------- ---------- ------ -------- --------- ------- ------- ----------------------------------------------------------------------------- ------------------------ ---- ------- ------------------- -------
    1    2011-09-20 01:00:00  2011-09-20 01:00:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_00_00_34z_image_lev1.fits aia__lev1:171:1095555635  N/A      No 2017-08-03 19:41:00 66200.0
    2    2011-09-20 01:00:00  2011-09-20 01:00:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_00_00_34z_image_lev1.fits aia__lev1:171:1095555635  N/A      No 2017-08-03 19:41:00 66200.0
    3    2011-09-20 01:15:00  2011-09-20 01:15:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_15_00_34z_image_lev1.fits aia__lev1:171:1095556535  N/A      No 2017-08-03 19:41:00 66200.0
    4    2011-09-20 01:15:00  2011-09-20 01:15:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_15_00_34z_image_lev1.fits aia__lev1:171:1095556535  N/A      No 2017-08-03 19:41:00 66200.0
    5    2011-09-20 01:30:00  2011-09-20 01:30:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_30_00_34z_image_lev1.fits aia__lev1:171:1095557435  N/A      No 2017-08-03 19:41:01 66200.0
    6    2011-09-20 01:30:00  2011-09-20 01:30:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_30_00_34z_image_lev1.fits aia__lev1:171:1095557435  N/A      No 2017-08-03 19:41:01 66200.0
    7    2011-09-20 01:45:00  2011-09-20 01:45:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_45_00_34z_image_lev1.fits aia__lev1:171:1095558335  N/A      No 2017-08-03 19:41:01 66200.0
    8    2011-09-20 01:45:00  2011-09-20 01:45:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_45_00_34z_image_lev1.fits aia__lev1:171:1095558335  N/A      No 2017-08-03 19:41:01 66200.0
    9    2011-09-20 02:00:00  2011-09-20 02:00:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t02_00_00_34z_image_lev1.fits aia__lev1:171:1095559235  N/A      No 2017-08-03 19:41:01 66200.0
   10    2011-09-20 02:00:00  2011-09-20 02:00:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t02_00_00_34z_image_lev1.fits aia__lev1:171:1095559235  N/A      No 2017-08-03 19:41:01 66200.0


If you then do a second query::

  >>> db.fetch(a.Time("2011-09-20T01:00:00", "2011-09-20T02:15:00"),
               a.Instrument('AIA'), a.vso.Sample(15*u.min))
  >>> db.commit()
  >>> db
  <Table length=12>
  id  observation_time_start observation_time_end instrument source provider  physobs  wavemin wavemax                                      path                                              fileid          tags starred    download_time      size
  str2         str19                 str19            str3     str3    str4      str9     str4    str4                                      str77                                              str24           str3   str2         str19          str7
  ---- ---------------------- -------------------- ---------- ------ -------- --------- ------- ------- ----------------------------------------------------------------------------- ------------------------ ---- ------- ------------------- -------
    1    2011-09-20 01:00:00  2011-09-20 01:00:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_00_00_34z_image_lev1.fits aia__lev1:171:1095555635  N/A      No 2017-08-03 19:41:00 66200.0
    2    2011-09-20 01:00:00  2011-09-20 01:00:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_00_00_34z_image_lev1.fits aia__lev1:171:1095555635  N/A      No 2017-08-03 19:41:00 66200.0
    3    2011-09-20 01:15:00  2011-09-20 01:15:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_15_00_34z_image_lev1.fits aia__lev1:171:1095556535  N/A      No 2017-08-03 19:41:00 66200.0
    4    2011-09-20 01:15:00  2011-09-20 01:15:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_15_00_34z_image_lev1.fits aia__lev1:171:1095556535  N/A      No 2017-08-03 19:41:00 66200.0
    5    2011-09-20 01:30:00  2011-09-20 01:30:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_30_00_34z_image_lev1.fits aia__lev1:171:1095557435  N/A      No 2017-08-03 19:41:01 66200.0
    6    2011-09-20 01:30:00  2011-09-20 01:30:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_30_00_34z_image_lev1.fits aia__lev1:171:1095557435  N/A      No 2017-08-03 19:41:01 66200.0
    7    2011-09-20 01:45:00  2011-09-20 01:45:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_45_00_34z_image_lev1.fits aia__lev1:171:1095558335  N/A      No 2017-08-03 19:41:01 66200.0
    8    2011-09-20 01:45:00  2011-09-20 01:45:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t01_45_00_34z_image_lev1.fits aia__lev1:171:1095558335  N/A      No 2017-08-03 19:41:01 66200.0
    9    2011-09-20 02:00:00  2011-09-20 02:00:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t02_00_00_34z_image_lev1.fits aia__lev1:171:1095559235  N/A      No 2017-08-03 19:41:01 66200.0
   10    2011-09-20 02:00:00  2011-09-20 02:00:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t02_00_00_34z_image_lev1.fits aia__lev1:171:1095559235  N/A      No 2017-08-03 19:41:01 66200.0
   11    2011-09-20 02:15:00  2011-09-20 02:15:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t02_15_00_34z_image_lev1.fits aia__lev1:171:1095560135  N/A      No 2017-08-03 19:42:19 66200.0
   12    2011-09-20 02:15:00  2011-09-20 02:15:01        AIA    SDO     JSOC intensity    17.1    17.1 /home/stuart/sunpy/data/aia_lev1_171a_2011_09_20t02_15_00_34z_image_lev1.fits aia__lev1:171:1095560135  N/A      No 2017-08-03 19:42:19 66200.0


A query can then be performed against the database to get the records.

  >>> entries = db.query(a.Time("2011-09-20T01:45:00", "2011-09-20T02:15:00"), a.Instrument('AIA'))
  >>> len(entries)
  4

You can see that only two extra records were added to the database. For more
information check out the :ref:`database_guide`.
