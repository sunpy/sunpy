.. _sunpy-coordinates:

SunPy coordinates
*****************

This sub-package contains:

* A robust framework for working with solar-physics coordinate systems
* Functions to obtain the locations of solar-system bodies (`sunpy.coordinates.ephemeris`)
* Functions to calculate Sun-specific coordinate information (`sunpy.coordinates.sun`)

The SunPy coordinate framework extends the
:ref:`Astropy coordinates framework <astropy:astropy-coordinates>`.

Supported Coordinate Systems
============================

.. list-table::
   :widths: auto
   :header-rows: 1

   * - Coordinate system
     - Abbreviation
     - SunPy/Astropy equivalent
     - Notes
   * - Heliocentric Aries Ecliptic (Mean)
     - HAE (also HEC)
     - Astropy's `~astropy.coordinates.HeliocentricMeanEcliptic`
     -
   * - Heliocentric Cartesian
     - HCC
     - `~sunpy.coordinates.frames.Heliocentric`
     -
   * - Heliocentric Earth Ecliptic
     - HEE
     - `~sunpy.coordinates.frames.HeliocentricEarthEcliptic`
     -
   * - Heliocentric Earth Equatorial
     - HEEQ (also HEQ)
     - `~sunpy.coordinates.frames.HeliographicStonyhurst`
     - Use a Cartesian representation
   * - Heliocentric Inertial
     - HCI
     - `~sunpy.coordinates.frames.HeliocentricInertial`
     -
   * - Heliocentric Radial
     - HCR
     - similar to `~sunpy.coordinates.frames.Heliocentric`
     - Use a cylindrical representation, *but* with a 90-degree offset in ``psi``
   * - Heliocentric/Heliographic Radial-Tangential-Normal
     - HGRTN
     - similar to `~sunpy.coordinates.frames.Heliocentric`
     - The axes are permuted, with HCC X, Y, Z equivalent respectively to HGRTN Y, Z, X
   * - Heliographic Carrington
     - HGC
     - `~sunpy.coordinates.frames.HeliographicCarrington`
     -
   * - Heliographic Stonyhurst
     - HGS
     - `~sunpy.coordinates.frames.HeliographicStonyhurst`
     -
   * - Helioprojective Cartesian
     - HPC
     - `~sunpy.coordinates.frames.Helioprojective`
     -
   * - Geocentric Earth Equatorial (Mean)
     - GEI
     - `~sunpy.coordinates.frames.GeocentricEarthEquatorial`
     -
   * - Geographic
     - GEO
     - Astropy's `~astropy.coordinates.ITRS`
     - The precise geographic definitions may differ
   * - Geocentric Solar Ecliptic
     - GSE
     - `~sunpy.coordinates.frames.GeocentricSolarEcliptic`
     -


For a description of these coordinate systems,
see `Thompson (2006) <https://doi.org/10.1051/0004-6361:20054262>`_
and `Franz & Harper (2002) <https://doi.org/10.1016/S0032-0633(01)00119-2>`_
(and `corrected version <https://www2.mps.mpg.de/homes/fraenz/systems/systems3art/systems3art.html>`_).


Getting Started
===============

The easiest interface to work with coordinates is through the `~astropy.coordinates.SkyCoord` class::

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord
  >>> from sunpy.coordinates import frames

  >>> c = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame=frames.Helioprojective)
  >>> c = SkyCoord(x=-72241.0*u.km, y=361206.1*u.km, z=589951.4*u.km, frame=frames.Heliocentric)
  >>> c = SkyCoord(70*u.deg, -30*u.deg, frame=frames.HeliographicStonyhurst)
  >>> c
  <SkyCoord (HeliographicStonyhurst: obstime=None): (lon, lat, radius) in (deg, deg, km)
      (70., -30., 695700.)>


It is also possible to use strings to specify the frame but in that case make sure to
explicitly import `sunpy.coordinates` as it registers SunPy's coordinate frames with
the Astropy coordinates framework::

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord

  >>> import sunpy.coordinates
  >>> c = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame='helioprojective', observer='earth')
  >>> c
  <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=earth): (Tx, Ty) in arcsec
      (-100., 500.)>


`~astropy.coordinates.SkyCoord` and all coordinate frames
support array coordinates. These work the same as single-value coordinates,
but they store multiple coordinates in a single object. When you're going to
apply the same operation to many different coordinates, this is a better choice
than a list of `~astropy.coordinates.SkyCoord` objects, because it will be
*much* faster than applying the operation to each
`~astropy.coordinates.SkyCoord` in a ``for`` loop::

   >>> c = SkyCoord([-500, 400]*u.arcsec, [100, 200]*u.arcsec, frame=frames.Helioprojective)
   >>> c
   <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=None): (Tx, Ty) in arcsec
       [(-500.,  100.), ( 400.,  200.)]>
   >>> c[0]
   <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=None): (Tx, Ty) in arcsec
       (-500.,  100.)>


Accessing Coordinates
---------------------

Individual coordinates can be accessed via attributes on the SkyCoord object,
but the names of the components of the coordinates can depend on the the frame and the chosen
representation (e.g., Cartesian versus spherical).

`~sunpy.coordinates.Helioprojective`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the helioprojective frame, the theta_x and theta_y components are accessed as
``Tx`` and ``Ty``, respectively::

  >>> c = SkyCoord(-500*u.arcsec, 100*u.arcsec, frame=frames.Helioprojective)
  >>> c.Tx
  <Longitude -500. arcsec>
  >>> c.Ty
  <Latitude 100. arcsec>

`~sunpy.coordinates.Heliocentric`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Heliocentric is typically used with Cartesian components::

  >>> c = SkyCoord(-72241.0*u.km, 361206.1*u.km, 589951.4*u.km, frame=frames.Heliocentric)
  >>> c.x
  <Quantity -72241. km>
  >>> c.y
  <Quantity 361206.1 km>
  >>> c.z
  <Quantity 589951.4 km>

`~sunpy.coordinates.HeliographicStonyhurst` and `~sunpy.coordinates.HeliographicCarrington`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both of the heliographic frames have the components of latitude, longitude and radius::

   >>> c = SkyCoord(70*u.deg, -30*u.deg, frame=frames.HeliographicStonyhurst)
   >>> c.lat
   <Latitude -30. deg>
   >>> c.lon
   <Longitude 70. deg>
   >>> c.radius
   <Distance 695700. km>

Heliographic Stonyhurst, when used with Cartesian components, is known as Heliocentric
Earth Equatorial (HEEQ).  Here's an example of how to use
`~sunpy.coordinates.frames.HeliographicStonyhurst` for HEEQ coordinates::

  >>> c = SkyCoord(-72241.0*u.km, 361206.1*u.km, 589951.4*u.km,
  ...              representation_type='cartesian', frame=frames.HeliographicStonyhurst)
  >>> c.x
  <Quantity -72241. km>
  >>> c.y
  <Quantity 361206.1 km>
  >>> c.z
  <Quantity 589951.4 km>

Transforming Between Coordinate Frames
======================================

Both `~astropy.coordinates.SkyCoord` and
`~astropy.coordinates.BaseCoordinateFrame` instances have a
`~astropy.coordinates.SkyCoord.transform_to` method. This can be used to
transform the frame to any other frame, either implemented in SunPy or in
Astropy (see also :ref:`astropy-coordinates-transforming`).
An example of transforming the center of the solar disk to Carrington
coordinates is::

   >>> c = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective, obstime="2017-07-26",
   ...              observer="earth")
   >>> c
   <SkyCoord (Helioprojective: obstime=2017-07-26T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty) in arcsec
       (0., 0.)>
   >>> c.transform_to(frames.HeliographicCarrington)
   <SkyCoord (HeliographicCarrington: obstime=2017-07-26T00:00:00.000, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (lon, lat, radius) in (deg, deg, AU)
       (283.95956776, 5.31701821, 0.00465047)>

It is also possible to transform to any coordinate system implemented in Astropy. This can be used to find the position of the solar limb in AltAz equatorial coordinates::

    >>> from astropy.coordinates import EarthLocation, AltAz
    >>> time = '2017-07-11 15:00'
    >>> greenbelt = EarthLocation(lat=39.0044*u.deg, lon=-76.8758*u.deg)
    >>> greenbelt_frame = AltAz(obstime=time, location=greenbelt)
    >>> west_limb = SkyCoord(900*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective,
    ...                      observer=greenbelt.get_itrs(greenbelt_frame.obstime), obstime=time)  # doctest: +REMOTE_DATA
    >>> west_limb.transform_to(greenbelt_frame)  # doctest: +REMOTE_DATA
    <SkyCoord (AltAz: obstime=2017-07-11 15:00:00.000, location=(1126916.53031967, -4833386.58391627, 3992696.62211575) m, pressure=0.0 hPa, temperature=0.0 deg_C, relative_humidity=0.0, obswl=1.0 micron): (az, alt, distance) in (deg, deg, m)
        (111.40782056, 57.1660434, 1.51859559e+11)>

Observer Location Information
=============================

The `~sunpy.coordinates.frames.Helioprojective`, `~sunpy.coordinates.frames.Heliocentric`
and `~sunpy.coordinates.frames.HeliographicCarrington` frames are defined by the location of
the observer. For example in `~sunpy.coordinates.frames.Helioprojective` the
observer is at the origin of the coordinate system. This information is encoded
in such observer-based frames as the ``observer`` attribute, which is itself an instance of the
`~sunpy.coordinates.frames.HeliographicStonyhurst` frame.
The ``observer`` can be a string for a solar-system body (e.g., "earth" or
"mars"), and
`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst` will be used
with the specified ``obstime`` to fully specify the observer location.  If
the observer location is not fully specified, or not present at all, most
transformations cannot be performed.
The location of the observer is automatically populated from meta data when
coordinate frames are created using map.

In the case of `~sunpy.coordinates.frames.HeliographicCarrington`, one can specify ``observer='self'`` to indicate that the coordinate itself should be used as the observer for defining the coordinate frame.

It is possible to convert from a `~sunpy.coordinates.frames.Helioprojective`
frame with one observer location to another
`~sunpy.coordinates.frames.Helioprojective` frame with a different observer
location, by converting through `~sunpy.coordinates.frames.Heliographic`, this
does involve making an assumption of the radius of the Sun to calculate the
position on the solar sphere. The conversion can be performed as follows::

  # Input coordinate
  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, observer="earth", obstime="2017-07-26", frame=frames.Helioprojective)

  Define a new Helioprojective frame with a different observer.
  >>> import sunpy.coordinates
  >>> hpc_out = sunpy.coordinates.Helioprojective(observer="venus", obstime="2017-07-26")

  Perform the transformation from one to the other.
  >>> hpc2 = hpc1.transform_to(hpc_out)

An example with two maps, named ``aia`` and ``stereo``::

  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=aia.coordinate_frame)  # doctest: +SKIP
  >>> hpc2 = hpc1.transform_to(stereo.coordinate_frame)  # doctest: +SKIP


Design of the Coordinates Sub-Package
=====================================

This sub-package works by defining a collection of ``Frames``
(`sunpy.coordinates.frames`), which exists on a transformation graph, where the
transformations between the coordinate frames are then defined and registered
with the transformation graph (`sunpy.coordinates.transformations`). It is also
possible to transform SunPy frames to Astropy frames.

Positions within these ``Frames`` are stored as a ``Representation`` of a
coordinate, a representation being a description of a point in a Cartesian,
spherical or cylindrical system (see :ref:`astropy-coordinates-representations`). A frame
that contains a representation of one or many points is said to have been
'realized'.

For a more in depth look at the design and concepts of the Astropy coordinates
system see :ref:`astropy-coordinates-overview`



Frames and SkyCoord
-------------------

The `~astropy.coordinates.SkyCoord` class is a high level wrapper around the
`astropy.coordinates` sub-package. It provides an easier way to create and transform
coordinates, by using string representations for frames rather than the classes
themselves and some other usability improvements, for more information see the
`~astropy.coordinates.SkyCoord` documentation.

The main advantage provided by `~astropy.coordinates.SkyCoord` is the support it
provides for caching Frame attributes. Frame attributes are extra data specified
with a frame, some examples in `sunpy.coordinates` are ``obstime`` or
``observer`` for observer location. Only the frames where this data is
meaningful have these attributes, i.e. only the Helioprojective frames have
``observer``. However, when you transform into another frame and then back to a
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

More Detailed Information
=========================

.. toctree::
   :maxdepth: 1

   carrington
   rotatedsunframe
   velocities
   wcs
   other_api


Reference/API
=============

.. automodapi:: sunpy.coordinates

.. automodapi:: sunpy.coordinates.ephemeris

.. automodapi:: sunpy.coordinates.sun

.. automodapi:: sunpy.coordinates.utils
    :no-inheritance-diagram:


Attribution
===========

Some of this documentation was adapted from Astropy under the terms of the `BSD
License
<https://raw.githubusercontent.com/astropy/astropy/master/LICENSE.rst>`_.

This sub-package was initially developed by Pritish Chakraborty as part of GSOC 2014 and Stuart Mumford.
