SunPy coordinates
=================

The SunPy coordinates submodule is an implementation of the common solar physics
coordinate frames using the :ref:`Astropy coordinates framework <astropy:astropy-coordinates>`.


Getting Started
---------------

The easiest interface to the coordinates module is through the `~astropy.coordinates.SkyCoord` class::

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord
  >>> from sunpy.coordinates import frames

  >>> c = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame=frames.Helioprojective)
  >>> c = SkyCoord(x=-72241.0*u.km, y=361206.1*u.km, z=589951.4*u.km, frame=frames.Heliocentric)
  >>> c = SkyCoord(70*u.deg, -30*u.deg, frame=frames.HeliographicStonyhurst)
  >>> c
  <SkyCoord (HeliographicStonyhurst: obstime=None): (lon, lat, rad) in (deg, deg, km)
      (70.0, -30.0, 695508.0)>


It is also possible to use strings to define the frame but in that case make sure to
explicitly import `sunpy.coordinates` as it registers solar coordinate frames with
astropy coordinates.::

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord

  >>> import sunpy.coordinates
  >>> c = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame='helioprojective')
  >>> c
  <SkyCoord (Helioprojective: D0=149597870.7 km, obstime=None, L0=0.0 deg, B0=0.0 deg, rsun=695508.0 km): (Tx, Ty) in arcsec
    (-100.,  500.)>


SunPy implements support for the following solar physics coordinate systems:

* Helioprojective (Cartesian) `~sunpy.coordinates.frames.Helioprojective`
* Heliocentric `~sunpy.coordinates.frames.Heliocentric`
* Heliographic Stonyhurst `~sunpy.coordinates.frames.HeliographicStonyhurst`
* Heliographic Carrington `~sunpy.coordinates.frames.HeliographicCarrington`

for a complete description of these frames see `sunpy.coordinates.frames`, for
a more detailed description of the frames see `Thompson (2006) <http://dx.doi.org/10.1051/0004-6361:20054262>`_


`~astropy.coordinates.SkyCoord` and all other `~astropy.coordinates` objects
also support array coordinates. These work the same as single-value coordinates,
but they store multiple coordinates in a single object. When you're going to
apply the same operation to many different coordinates, this is a better choice
than a list of `~astropy.coordinates.SkyCoord` objects, because it will be
*much* faster than applying the operation to each
`~astropy.coordinates.SkyCoord` in a ``for`` loop.
::

   >>> c = SkyCoord([-500, 400]*u.arcsec, [100, 200]*u.arcsec, frame=frames.Helioprojective)
   >>> c
   <SkyCoord (Helioprojective: obstime=None, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=None): (lon, lat, radius) in (deg, deg, AU)
       ( 0.,  0.,  1.)>): (Tx, Ty) in arcsec
       [(-500.,  100.), ( 400.,  200.)]>
   >>> c[0]
   <SkyCoord (Helioprojective: obstime=None, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=None): (lon, lat, radius) in (deg, deg, AU)
       ( 0.,  0.,  1.)>): (Tx, Ty) in arcsec
       (-500.,  100.)>


Accessing Coordinates
^^^^^^^^^^^^^^^^^^^^^

Individual coordinates can be accessed via attributes on the SkyCoord object,
but the names of the components of the coordinates for each frame differ. For a full
description of all the properties of the frames see `sunpy.coordinates.frames`.

``HelioProjective``
###################

For the helioprojective frame the coordinates are access as ``Tx`` and ``Ty`` representing theta x and y. These are the same coordinates that are often referred to as 'solar-x' and 'solar-y'.::

  >>> c = SkyCoord(-500*u.arcsec, 100*u.arcsec, frame=frames.Helioprojective)
  >>> c.Tx
  <Longitude180 -500.0 arcsec>
  >>> c.Ty
  <Latitude 100.0 arcsec>

``Heliocentric``
################

Heliocentric normally a Cartesian frame so the coordinates are accessed as ``x,y,z``::

  >>> c = SkyCoord(-72241.0*u.km, 361206.1*u.km, 589951.4*u.km, frame=frames.Heliocentric)
  >>> c.x
  <Quantity -72241.0 km>
  >>> c.y
  <Quantity 361206.1 km>
  >>> c.z
  <Quantity 589951.4 km>

``HeliographicStonyhurst`` and ``HeliographicCarrington``
#########################################################

Both the heliographic frames use latitude, longitude and radius which are accessed as follows::

   >>> c = SkyCoord(70*u.deg, -30*u.deg, frame=frames.HeliographicStonyhurst)
   >>> c.lat
   <Latitude -30.0 deg>
   >>> c.lon
   <Longitude180 70.0 deg>
   >>> c.radius
   <Distance 695508.0 km>

Transforming Between Coordinate Frames
--------------------------------------

Both `~astropy.coordinates.SkyCoord` and `~astropy.coordinates.BaseCoordinateFrame` instances have a `~astropy.coordinates.SkyCoord.transform_to` method. This can be used to transform the frame to any other frame, either implemented in SunPy or in Astropy. An example of transforming the center of the solar disk to Carrington coordinates is::

   >>> c = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective, obstime="2017-07-26")
   >>> c
   <SkyCoord (Helioprojective: obstime=2017-07-26 00:00:00, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2017-07-26 00:00:00): (lon, lat, radius) in (deg, deg, AU)
       ( 0.,  5.31701821,  1.01567428)>): (Tx, Ty) in arcsec
       ( 0.,  0.)>

   >>> c.transform_to(frames.HeliographicCarrington)
   <SkyCoord (HeliographicCarrington: obstime=2017-07-26 00:00:00): (lon, lat, radius) in (deg, deg, km)
       (-76.00701638,  5.31701821,  695508.00000058)>

It is also possible to transform to any coordinate system implemented in Astropy. This can be used to find the position of the solar limb in AltAz equatorial coordinates::

    >>> from astropy.coordinates import EarthLocation, AltAz

    >>> time = '2017-07-11 15:00'
    >>> greenbelt = EarthLocation(lat=39.0044*u.deg, lon=-76.8758*u.deg)
    >>> greenbelt_frame = AltAz(obstime=time, location=greenbelt)

    >>> west_limb = SkyCoord(900*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective, obstime=time)
    >>> west_limb.transform_to(greenbelt_frame)
    <AltAz Coordinate (obstime=2017-07-11 15:00:00.000, location=(1126916.53031967, -4833386.58391627, 3992696.622115747) m, pressure=0.0 hPa, temperature=0.0 deg_C, relative_humidity=0, obswl=1.0 micron): (az, alt, distance) in (deg, deg, m)
        ( 111.40839171,  57.16645763,   1.51860261e+11)>


Observer Location Information
-----------------------------

Both `~sunpy.coordinates.frames.Helioprojective` and
`~sunpy.coordinates.frames.Heliocentric` frames are defined by the location of
the observer. For example in `~sunpy.coordinates.frames.Helioprojective` the
observer is at the origin of the coordinate system. This information is encoded
in the `~sunpy.coordinates.frames.Helioprojective` and
`~sunpy.coordinates.frames.Heliocentric` frames as the ``observer`` attribute,
which is itself an instance of the
`~sunpy.coordinates.frames.HeliographicStonyhurst` frame. The default observer
location is set to the position of the Earth (using
`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`) as long as the
``obstime`` attribute is specified. If the ``obstime`` attribute is not set the
observer defaults to ``(0°, 0°, 1 AU)`` i.e. the mean position of the Earth. The
location of the observer is automatically populated from meta data when
coordinate frames are created using map.

It is possible to convert from a `~sunpy.coordinates.frames.Helioprojective`
frame with one observer location to another
`~sunpy.coordinates.frames.Helioprojective` frame with a different observer
location, by converting through `~sunpy.coordinates.frames.Heliographic`, this
does involve making an assumption of the radius of the Sun to calculate the
position on the solar sphere. The conversion can be performed as follows::

  # Input coordinate
  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, observer="earth", obstime="2017-07-26", frame=frames.Helioprojective)
  # Define a new Helioprojective frame with a different observer.
  >>> hpc_out = sunpy.coordinates.Helioprojective(observer="venus", obstime="2017-07-26")
  # Perform the transformation from one to the other.
  >>> hpc2 = hpc1.transform_to(hpc_out)

An example with two maps, i.e. ``aia`` and ``stereo``::

  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=aia.coordinate_frame)
  >>> hpc2 = hpc1.transform_to(stereo.coordinate_frame)


Design of the Coordinates Module
--------------------------------

This module works by defining a collection of ``Frames``
(`sunpy.coordinates.frames`), which exists on a transformation graph, where the
transformations between the coordinate frames are then defined and registered
with the transformation graph (`sunpy.coordinates.transformations`). It is also
possible to transform SunPy frames to Astropy frames.

Positions within these ``Frames`` are stored as a ``Representation`` of a
coordinate, a representation being a description of a point in a Cartesian,
spherical or cylindrical system (`sunpy.coordinates.representation`). A frame
that contains a representation of one or many points is said to have been
'realized'.

For a more in depth look at the design and concepts of the Astropy coordinates
system see :ref:`astropy-coordinates-overview`



Frames and SkyCoord
^^^^^^^^^^^^^^^^^^^

The `~astropy.coordinates.SkyCoord` class is a high level wrapper around the
`astropy.coordinates` package. It provides an easier way to create and transform
coordinates, by using string representations for frames rather than the classes
themselves and some other usability improvements, for more information see the
`~astropy.coordinates.SkyCoord` documentation.

The main advantage provided by `~astropy.coordinates.SkyCoord` is the support it
provides for caching Frame attributes. Frame attributes are extra data specified
with a frame, some examples in `sunpy.coordinates` are ``obstime`` or
``observer`` for observer location. Only the frames where this data is
meaningful have these attributes, i.e. only the Helioprojective frames have
``observe``. However, when you transform into another frame and then back to a
projective frame using `SkyCoord` it will remember the attributes previously
provided, and repopulate the final frame with them. If you were to do
transformations using the Frames alone this would not happen.

The most important implication for this in `sunpy.coordinates` is the ``rsun``
parameter in the projective frames. If you create a projective frame with a
``rsun`` attribute, if you convert back to a projective frame it will be set
correctly. It should also be noted that, if you create a Heliographic frame and
then transform to a projective frame with an ``rsun`` attribute, it will not
match the ``radius`` coordinate in the Heliographic frame. This is because you may
mean to be describing a point above the defined 'surface' of the Sun.


Coordinates and WCS
-------------------

The `sunpy.coordinates` package provides a mapping between FITS-WCS CTYPE
convention and the coordinate frames as defined in `sunpy.coordinates`. This is
used via the `astropy.wcs.utils.wcs_to_celestial_frame` function, with which the
SunPy frames are registered upon being imported. This list is used by packages
such as ``wcsaxes`` to convert from `astropy.wcs.WCS` objects to coordinate
frames.

The `sunpy.map.GenricMap` class creates `astropy.wcs.WCS` objects as
``amap.wcs``, however, it adds some extra attributes to the `~astropy.wcs.WCS`
object to be able to fully specify the coordinate frame. It adds
``heliographic_observer`` and ``rsun``.

If you want to obtain a un-realized coordinate frame corresponding to a
`~sunpy.map.GenericMap` object you can do the following::

  >>> import sunpy.map
  >>> from sunpy.data.sample import AIA_171_IMAGE

  >>> amap = sunpy.map.Map(AIA_171_IMAGE)
  >>> amap.observer_coordinate
  <Helioprojective Frame (obstime=2011-06-07 06:33:02.770000, rsun=696000000.0 m, observer=<HeliographicStonyhurst Coordinate (obstime=None): (lon, lat, radius) in (deg, deg, m)
      ( 0.,  0.048591,   1.51846026e+11)>)>


which is equivalent to::

  >>> from astropy.wcs.utils import wcs_to_celestial_frame
  >>> wcs_to_celestial_frame(amap.wcs)
  <Helioprojective Frame (obstime=2011-06-07 06:33:02.770000, rsun=696000000.0 m, observer=<HeliographicStonyhurst Coordinate (obstime=None): (lon, lat, radius) in (deg, deg, m)
      ( 0.,  0.048591,   1.51846026e+11)>)>


.. automodapi:: sunpy.coordinates
    :headings: ^#

.. automodapi:: sunpy.coordinates.frames
    :headings: ^#

.. automodapi:: sunpy.coordinates.transformations
    :headings: ^#

.. automodapi:: sunpy.coordinates.ephemeris
    :headings: ^#

.. automodapi:: sunpy.coordinates.representation
    :headings: ^#

.. automodapi:: sunpy.coordinates.offset_frame
    :headings: ^#

.. automodapi:: sunpy.coordinates.wcs_utils
    :headings: ^#


Attribution
-----------

Some of this documentation was adapted from Astropy under the terms of the `BSD
License
<https://raw.githubusercontent.com/astropy/astropy/master/licenses/LICENSE.rst>`_.

This package was initially developed by Pritish Chakraborty as part of GSOC 2014 and Stuart Mumford.
