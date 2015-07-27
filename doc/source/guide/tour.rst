---------------------
A brief tour of SunPy
---------------------

This brief tutorial will walk you through some
of the functionality offered by SunPy. Start by reading this tutorial
and trying out some of the examples demonstrated. Once you've completed the
tutorial check out the rest of the :doc:`User Guide </guide/index>` for a more
thorough look at the functionality available.

Sample Data
-----------
This tour makes use of a number of sample data files which you will need to
download. If have not already done so please follow the instruction here :ref:`sample-data`.

Maps
----
Maps are the primary data type in SunPy they are spatially and / or temporally aware
data arrays. There are maps for a 2D image, a time series of 2D images or temporally aligned 2D images.

**Creating a Map**

SunPy supports many different data products from various sources 'out of the box'. We
shall use SDO's AIA instrument as an example in this tutorial. The general way to create
a Map from one of the supported data products is with the `~sunpy.map.map()` function from the `~sunpy.map` submodule.
`~sunpy.map.map()` takes either a filename, a list of filenames or a data array and header. We can test map with:

.. plot::
    :include-source:

    import sunpy.data.sample
    import sunpy.map
    aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    aia.peek()

This returns a map named `aia` which can be manipulated with standard SunPy map commands.
For more information about maps checkout the :doc:`map guide <data_types/maps>`
and the :ref:`map`.

Lightcurve
----------

SunPy handles time series data, fundamental to the study of any real world phenomenon,
by creating a lightcurve object. A lightcurve consists of two parts; times and measurements taken at those times. The
data can either be in your current Python session, alternatively within a local or
remote file. Let's create some fake data and pass it into a lightcurve object.

.. plot::
    :include-source:

    import sunpy.data.sample
    from sunpy.lightcurve import LightCurve
    times = np.arange(1000) * 2.0
    signal = np.sin(np.arange(1000)*0.02 ) + np.random.random(1000)
    light_curve = LightCurve.create({"signal": signal},index = times)
    light_curve.peek()

Within LightCurve.create, we have a dictionary that contains a single entry with key
"signal" containing a list of 1000 entries (0-999). The accompanying set of times is
passed in via the index keyword argument. If no times are passed into index, a default
set of time indices is generated.

For more information about lightcurves, check out the
:doc:`lightcurve guide <data_types/lightcurve>` and the
and the :ref:`lightcurve_code_ref`.

.. this should be a better example, for example grabbing goes data...

Spectra
-------

SunPy has spectral support for instruments which have such a capacity. CALLISTO,
an international network of Solar Radio Spectrometers, is a specific example.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import sunpy.spectra
    import sunpy.data.sample
    from sunpy.spectra.sources.callisto import CallistoSpectrogram
    image = CallistoSpectrogram.read(sunpy.data.sample.CALLISTO_IMAGE)
    image.peek()

For more information about spectra, check out the :doc:`spectra guide <data_types/spectra>`
and the :ref:`spectra_code_ref`.

Plotting
--------

SunPy uses a matplotlib like interface to it's plotting so more complex
plots can be built by combining SunPy with matplotlib.

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
a red colormap, a colorbar on the right-hand side and a title and some
labels.

There is lot going on here, but we will walk you through the example. Briefly,
the first line is just importing SunPy. On the second line we create a
SunPy Map object which is basically just a spatially-aware image or data array.
On the last line we then plot the map object, using the built in 'quick plot' function `peek()`.

SunPy uses a matplotlib like interface to it's plotting so more complex
plots can be built by combining SunPy with matplotlib.

.. plot::
    :include-source:

    import sunpy.map
    import matplotlib.pyplot as plt
    import sunpy.data.sample
    aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    fig = plt.figure()
    ax = plt.subplot(111)
    aia.plot()
    aia.draw_limb()
    aia.draw_grid()
    plt.colorbar()
    aia.draw_limb()
    plt.show()

For more information check out :ref:`plotting`.

Solar Physical Constants
------------------------

SunPy contains a convenient list of solar-related physical constants. Here is
a short bit of code to get you started: ::

    >>> from sunpy.sun import constants as con

    # one astronomical unit (the average distance between the Sun and Earth)
    >>> print con.au
      Name   = Astronomical Unit
      Value  = 1.495978707e+11
      Error  = 0.0
      Units  = m
      Reference = IAU 2012 Resolution B2

    # the solar radius
    >>> print con.radius
      Name   = Solar radius
      Value  = 695508000.0
      Error  = 26000.0
      Units  = m
      Reference = Allen's Astrophysical Quantities 4th Ed.

Not all constants have a shortcut assigned to them (as above). The rest of the constants
are stored in a dictionary. The following code grabs the dictionary and gets all of the
keys.::

    >>> solar_constants = con.physical_constants
    >>> solar_constants.keys()   # doctest: +NORMALIZE_WHITESPACE
    ['solar flux unit', 'surface area', 'average density', 'radius', 'surface
    gravity', 'ellipticity', 'visual magnitude', 'center density', 'average
    angular size', 'absolute magnitude', 'sunspot cycle', 'effective
    temperature', 'aphelion distance', 'mean energy production', 'mass
    conversion rate', 'average intensity', 'volume', 'metallicity', 'moment of
    inertia', 'escape velocity', 'perihelion distance', 'GM', 'oblateness',
    'mean distance', 'age', 'mass', 'luminosity', 'center temperature']

You can also use the following function to print out a table of all of the values
available. ::

    >>> con.print_all()   # doctest: +NORMALIZE_WHITESPACE
    Name                                 Value            Units    Error
    -------------------------------------------------------------------------------------
    solar flux unit                      1e-22      W / (Hz m2)    0
    surface area                     6.087e+18               m2    0
    average density                       1409          kg / m3    0
    radius                         695508000.0                m    26000.0
    surface gravity                        274            m / s    0
    ellipticity                          5e-05                     0
    visual magnitude                    -26.75                     0
    center density                    162200.0          kg / m3    0
    average angular size                959.63           arcsec    0
    absolute magnitude                    4.83                     0
    sunspot cycle                         11.4               yr    0
    effective temperature               5778.0                K    0
    aphelion distance                1.521e+11                m    0
    mean energy production           0.0001937           J / kg    0
    mass conversion rate          4300000000.0           kg / s    0
    average intensity               20090000.0      W / (m2 sr)    0
    volume                          1.4122e+27               m3    0
    metallicity                         0.0122                     0.0
    moment of inertia                  5.7e+54          kg / m2    0
    escape velocity                   617700.0            m / s    0
    perihelion distance              1.471e+11                m    0
    GM                             132712000.0         km3 / s2    0
    oblateness                            8.01          marcsec    0.14
    mean distance              1.495978707e+11                m    0.0
    age                           4600000000.0               yr    100000000.0
    mass                            1.9891e+30               kg    5e+25
    luminosity                       3.846e+26                W    5e+22
    center temperature              15710000.0                K    0

These constants are provided as a convenience so that everyone is using the same
(accepted values). For more information check out :ref:`sun_code_ref`.

Quantities and Units
--------------------

Many capabilities in SunPy make use of physical quantities that are specified
with units. SunPy uses `~astropy.units` to
implement this functionality. For example, the solar radius above is a physical quantity
that can be expressed in length units.  In the example above ::

    from sunpy.sun import constants as con
    con.radius
    <Constant name=u'Solar radius' value=695508000.0 error=26000.0 units='m' reference=u"Allen's Astrophysical Quantities 4th Ed.">

shows the solar radius in units of meters.  It is simple to express the same physical quantity in different units::

    con.radius.to('km')
    <Quantity 695508.0 km>

To get the numerical value of the solar radius in kilometers - without the unit information - use ::

    con.radius.to('km').value
    695508.0

Quantities and units are simple and powerful tools for keeping track of the units you're working in, and make it
easy to convert the same physical quantity into different units.  To learn more about the capabilities of quantities
and units, please consult `the astropy tutorial <http://www.astropy.org/astropy-tutorials/Quantities.html>`__.
SunPy's approach to the adoption of quantities and units in the codebase is described
`here <https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0003.md>`__.

Here's a simple example of the power of units.  Suppose you have the radius of a circle and would like to calculate
its area.  The following code implements this ::

    >>> import numpy as np
    >>> import astropy.units as u
    >>> @u.quantity_input(radius=u.m)
    ... def circle_area(radius):
    ...     return np.pi * radius ** 2

The first line imports numpy, and the second line imports astropy's units module.  The beginning of the third line (the
"@" symbol) indicates that what follows is a Python decorator.  In this case, the decorator allows us to specify what
kind of unit the function input variable "radius" in the following function "circle_area" should have.  In this case,
it is meters.  The decorator checks that the input is convertible to the units specified in the decorator.  Calculating
the area of a circle with radius 4 meters using the function defined above is simple ::

    circle_area(4 * u.m)
    <Quantity 50.26548245743669 m2>

The units of the returned area are what we expect, namely the meters squared (m2).  However, we can also use other
units of measurement; for a circle with radius 4 kilometers ::

    circle_area(4 * u.km)
    <Quantity 50.26548245743669 km2>

Even although the input value of the radius was not in meters, the function does not crash; this is because the
input unit is convertible to meters.  This also works across different systems of measurement, for example ::

    circle_area(4 * u.imperial.foot)
    <Quantity 50.26548245743669 ft2>

However, if the input unit is not convertible to meters, then an error is thrown ::

    >>> circle_area(4 * u.second)   # doctest: +SKIP
    ...
    UnitsError: Argument 'radius' to function 'circle_area' must be in units convertable to 'm'.

Also, if no unit is specified, an error is thrown ::

    >>> circle_area(4)   # doctest: +SKIP
    ...
    TypeError: Argument 'radius' to function has 'circle_area' no 'unit' attribute. You may want to pass in an astropy Quantity instead.

Using units allows the user to be explicit about what the function
expects.  Units also make conversions very easy to do.  For example,
if you want the area of a circle in square feet, but were given
measurements in meters, then ::

    circle_area((4 * u.m).to(u.imperial.foot))
    <Quantity 541.0531502245425 ft2>

or ::

    >>> circle_area(4 * u.m).to(u.imperial.foot ** 2)
    <Quantity 541.0531502245426 ft2>

Astropy units and quantities are very powerful, and are used throughout SunPy.  To find out more about units and
quantities, please consult the `the astropy tutorial <http://www.astropy.org/astropy-tutorials/Quantities.html>`__ and
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


Getting at Data
---------------

Querying the VSO
----------------
There are a couple different ways to query and download data from the VSO using
SunPy. The method you should use depends first on your preference with respect
to query style: the main method of querying uses a syntax that is unique to
SunPy and may require some getting used to, but is extremely flexible and
powerful. A second
"legacy" API also exists which works is very much the same way as VSO_GET in
IDL.

Further, for each of the two query APIs there are interactive and
non-interactive versions available, depending on the type of work you are doing.

The below example demonstrates a simple query for SOHO EIT data using the
non-interactive version of the main API::

    >>> from sunpy.net import vso

    # create a new VSOClient instance
    >>> client = vso.VSOClient()

    # build our query
    >>> result = client.query(
    ...     vso.attrs.Time((2011, 9, 20, 1), (2011, 9, 20, 2)),
    ...     vso.attrs.Instrument('eit'))

    # print the number of matches
    >>> print("Number of records found: %d " % result.num_records())   # doctest: +NORMALIZE_WHITESPACE
    Number of records found: 4

    # download matches to /download/path
    >>> res = client.get(result, path="/download/path/{file}").wait()

Note that specifying a path is optional and if you do not specify one the files
will simply be downloaded into a temporary directory (e.g. /tmp/xyz).
For more information about vso client checkout the :doc:`vso guide <acquiring_data/vso>`.

Database Package
----------------

The database package offers the possibility to save retrieved data (e.g. via the
:mod:'sunpy.net.vso' package) onto a local or remote database. The database may be
a single file located on a local hard drive (if a SQLite database is used) or a
local or remote database server.
This makes it possible to fetch required data from the local database instead
of downloading it again from a remote server.

Querying a database is straightforward, as this example using VSO, shows. The example
demonstrates the useful feature which prevents storing the same data twice::


    >>> from sunpy.database import Database
    >>> from sunpy.net.vso.attrs import Time, Instrument
    >>> db = Database('sqlite:///')
    >>> entries = db.fetch(
    ...     Time('2012-08-05', '2012-08-05 00:00:05'),
    ...     Instrument('AIA'))
    >>> assert entries is None
    >>> len(db)
    4
    >>> entries = db.fetch(
    ...     Time('2012-08-05', '2012-08-05 00:00:05'),
    ...     Instrument('AIA'))
    >>> entries is None
    False
    >>> len(entries)
    4
    >>> len(db)
    4


Explanation: first, entries is None because the query has never been used for querying
the database -> query the VSO, add new entries to database, remember query hash.
In the second fetch, entries is not None because the query has already been used and
returns a list of database entries. For more information check out the :ref:`database_guide`.

Querying Helioviewer.org
------------------------

SunPy can be used to make several basic requests using the The `Helioviewer.org API <http://helioviewer.org/api/>`__
including generating a PNG and downloading a `JPEG 2000 <http://wiki.helioviewer.org/wiki/JPEG_2000>`__
image and loading it into a SunPy Map.


A simple example of a helioviewer query and a plot of the result follows.

.. plot::
    :include-source:

    from sunpy.net.helioviewer import HelioviewerClient
    import matplotlib.pyplot as plt
    from matplotlib.image import imread
    hv = HelioviewerClient()
    file = hv.download_png('2099/01/01', 4.8, "[SDO,AIA,AIA,304,1,100]", x0=0, y0=0, width=512, height=512)
    im = imread(file)
    plt.imshow(im)
    plt.axis('off')
    plt.show()

This downloads a PNG image of the latest AIA 304 image available on `Helioviewer.org <http://helioviewer.org>`_.  In the
 `~sunpy.net.helioviewer.HelioviewerClient.download_png` command the value, 4.8, refers to the image resolution in arcseconds per pixel (larger values mean lower resolution), x0 and y0 are the center points about which to focus and the width and height are the pixel values for the image dimensions. For more information checkout the :doc:`helioviewer guide <acquiring_data/helioviewer>`.
