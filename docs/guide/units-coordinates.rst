.. _units-coordinates-sunpy:

Units and Coordinates in SunPy
******************************

This section of the guide will talk about representing physical units and
physical coordinates in SunPy. SunPy makes use of :ref:`Astropy <astropy:astropy-coordinates>` for
both these tasks.


Units in SunPy
==============

All functions in SunPy that accept or return numbers associated with physcial
quantities accept and return `~astropy.units.Quantity` objects. These objects
represent a number (or an array of numbers) and a unit. This means SunPy is
always explicit about the units associated with a value. Quantities and units
are powerful tools for keeping track of variables with a physical meaning and
make it straightforward to convert the same physical quantity into different units.

In this section of the guide we will give a quick introduction to `astropy.units`
and then demostrate how to use units with SunPy.

To use units we must first import them from Astropy. To save on typing we usually
import units as ``u``::

   >>> import astropy.units as u

Once we have imported units we can create a quantity by multiplying a number by
a unit::

   >>> length = 10 * u.meter
   >>> length
   <Quantity 10. m>

A `~astropy.units.Quantity` has both a ``.unit`` and a ``.value`` attribute::

  >>> length.value
  10.0

  >>> length.unit
  Unit("m")

These `~astropy.units.Quantity` objects can also be converted to other units, or
unit systems::

  >>> length.to(u.km)
  <Quantity 0.01 km>

  >>> length.cgs
  <Quantity 1000. cm>

Probably most usefully, `~astropy.units.Quantity` objects will propogate units
through arithmetic operations when appropriate::

  >>> distance_start = 10 * u.mm
  >>> distance_end = 23 * u.km
  >>> length = distance_end - distance_start
  >>> length
  <Quantity 22.99999 km>

  >>> time = 15 * u.minute
  >>> speed = length / time
  >>> speed
  <Quantity 1.53333267 km / min>

However, operations which do not make physical sense for the units specified will cause an error::

  >>> length + time
  Traceback (most recent call last):
  ...
  astropy.units.core.UnitConversionError: Can only apply 'add' function to quantities with compatible dimensions


Quantities as function arguments
================================

An extremely useful addition to the base functionality of Quanitities is the ``@u.quantity_input`` decorator.
This allows you to specify required units for function arguments to ensure that the calculation within that
function always make physical sense. For instance, if we defined a function to calculate speed as above,
we might want the distance and time as inputs::

  >>> def speed(length, time):
  ...     return length / time

However, this requires that length and time both have the appropriate units. We therefore want to use
`~astropy.units.quantity_input` to enforce this, here we use
`function annotations <https://python-3-for-scientists.readthedocs.io/en/latest/python3_features.html#function-annotations>`__
to specify the units.

  >>> @u.quantity_input
  ... def speed(length: u.m, time: u.s):
  ...     return length / time

Now, when this function is called, if the units of length and time are not convertible to the units specified,
an error will be raised stating that the units are incorrect or missing::

  >>> speed(1*u.m, 10*u.m)
  Traceback (most recent call last):
  ...
  astropy.units.core.UnitsError: Argument 'time' to function 'speed' must be in units convertible to 's'.

  >>> speed(1*u.m, 10)
  ...
  Traceback (most recent call last):
  ...
  TypeError: Argument 'time' to function 'speed' has no 'unit' attribute. ... pass in an astropy Quantity instead.

Note that the units of the inputs do not have to be exactly the same as those in the function definition, as long
as they can be converted to those units. So for instance, passing in a time in minutes still works even though we
specified ``time: u.s``::

  >>> speed(1*u.m, 1*u.minute)
  <Quantity 1. m / min>

This may still not be quite as we want it, since we wanted the input time in seconds but the output is in m/min.
We can correct this by defining the function with an additional annotation::

  >>> @u.quantity_input
  ... def speed(length: u.m, time: u.s) -> u.m/u.s:
  ...     return length / time

This will force the output of the function to be converted to m/s before returning, so that you will always
have the same units on the output from this function::

  >>> speed(1*u.m, 1*u.minute)
  <Quantity 0.01666667 m / s>

Physical Coordinates in SunPy
=============================

In much the same way as `~astropy.units` are used for representing physical
quantities, SunPy uses `astropy.coordinates` to represent points in physical
space. This applies to both points in 3D space and projected coordinates in
images.

The astropy coordinates module is primarily used through the
`~astropy.coordinates.SkyCoord` class::

  >>> from astropy.coordinates import SkyCoord

To enable the use of the solar physics specific frames defined in SunPy we also
need to import them::

  >>> from sunpy.coordinates import frames

A SkyCoord object to represent a point on the Sun can then be created::

  >>> c = SkyCoord(70*u.deg, -30*u.deg, obstime="2017-08-01",
  ...              frame=frames.HeliographicStonyhurst)
  >>> c
  <SkyCoord (HeliographicStonyhurst: obstime=2017-08-01T00:00:00.000): (lon, lat, radius) in (deg, deg, km)
      (70., -30., 695700.)>

This `~astropy.coordinates.SkyCoord` object can then be transformed to any
other coordinate frame defined either in Astropy or SunPy, for example::

  >>> c.transform_to(frames.Helioprojective(observer="earth"))
  <SkyCoord (Helioprojective: obstime=2017-08-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, km)
      (769.96270814, -498.89715922, 1.51668773e+08)>


It is also possible to convert three dimensional positions to astrophysical
frames defined in Astropy, for example `~astropy.coordinates.ICRS`.

  >>> c.transform_to('icrs')
  <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, km)
    (49.84856512, 0.05394699, 1417743.94689472)>



Observer Location
-----------------

Both `~sunpy.coordinates.frames.Helioprojective` and
`~sunpy.coordinates.frames.Heliocentric` frames are defined based on the
position of the observer. Therefore to transform either of these frames to a
different frame the location of the observer must be known.
The observer can be specified for a coordinate object using the ``observer``
argument to `~astropy.coordinates.SkyCoord`.  For SunPy to calculate the
location of Earth or another solar-system body, it must know the time for
which the coordinate is valid; this is specified with the ``obstime`` argument.

Using the observer location it is possible to convert a coordinate as seen by
one observer to a coordinate seen by another::

  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, observer="earth",
  ...                 obstime="2017-07-26",
  ...                 frame=frames.Helioprojective)

  >>> hpc1.transform_to(frames.Helioprojective(observer="venus",
  ...                                          obstime="2017-07-26"))
  <SkyCoord (Helioprojective: obstime=2017-07-26T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'venus'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (-1285.47497992, 106.20918654, 0.72405937)>


Using Coordinates with SunPy Map
--------------------------------

.. plot::
   :include-source:

   SunPy Map uses coordinates to specify locations on the image, and to plot
   overlays on plots of maps. When a Map is created, a coordinate frame is
   constructed from the header information. This can be accessed using
   ``.coordinate_frame``:

   >>> import sunpy.map
   >>> from sunpy.data.sample import AIA_171_IMAGE   # doctest: +REMOTE_DATA
   >>> m = sunpy.map.Map(AIA_171_IMAGE)  # doctest: +REMOTE_DATA
   >>> m.coordinate_frame  # doctest: +REMOTE_DATA
   <Helioprojective Frame (obstime=2011-06-07T06:33:02.770, rsun=696000000.0 m, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770): (lon, lat, radius) in (deg, deg, m)
       (-0.00406308, 0.04787238, 1.51846026e+11)>)>

   This can be used when creating a `~astropy.coordinates.SkyCoord` object to set
   the coordinate system to that image:

   >>> from astropy.coordinates import SkyCoord
   >>> import astropy.units as u
   >>> c = SkyCoord(100 * u.arcsec, 10*u.arcsec, frame=m.coordinate_frame)  # doctest: +REMOTE_DATA
   >>> c  # doctest: +REMOTE_DATA
   <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.770, rsun=696000000.0 m, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770): (lon, lat, radius) in (deg, deg, m)
       (-0.00406308, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
       (100., 10.)>

   The `~astropy.coordinates.SkyCoord` object can be converted to a PixelPair object
   using `~sunpy.map.GenericMap.world_to_pixel`:

   >>> pixel_obj = m.world_to_pixel(c) # doctest: +REMOTE_DATA
   >>> pixel_obj # doctest: +REMOTE_DATA
   PixelPair(x=<Quantity 551.7680511 pix>, y=<Quantity 515.18266871 pix>)

   This `~astropy.coordinates.SkyCoord` object could also be used to plot a point
   on top of the map:

   >>> import matplotlib.pyplot as plt
   >>> ax = plt.subplot(projection=m)  # doctest: +SKIP
   >>> m.plot()  # doctest: +SKIP
   <matplotlib.image.AxesImage object at ...>
   >>> _ = ax.plot_coord(c, 'o')  # doctest: +SKIP

For more information on coordinates see :ref:`sunpy-coordinates` section of the :ref:`reference`.
