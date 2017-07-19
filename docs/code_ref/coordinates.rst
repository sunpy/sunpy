SunPy coordinates
=================

The SunPy coordinates submodule is an implementation of the common solar physics
coordinate frames using the :ref:`Astropy coordinates framework <astropy:astropy-coordinates>`.

.. note::

   The accuracy of the transformations in this module have not been rigorously
   verified. They have been compared to the previous `sunpy.wcs` implementation
   and match to numerical precision. Independent verification will be added at a
   later date.


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
astropy coordinates.

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord
  >>> import sunpy.coordinates
  >>> c = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame='helioprojective')
  >>> c
  <SkyCoord (Helioprojective: D0=149597870.7 km, obstime=None, L0=0.0 deg, B0=0.0 deg, rsun=695508.0 km): (Tx, Ty) in arcsec
    (-100.,  500.)>


SunPy implements support for the following solar physics coordinate systems:

* Helioprojective (Cartesian) `~sunpy.coordinates.frames.HelioProjective`
* Heliocentric `~sunpy.coordinates.frames.HelioCentric`
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
   <SkyCoord (HelioProjective: D0=149597870.7 km, obstime=None, L0=0.0 deg, B0=0.0 deg, rsun=695508.0 km): (Tx, Ty) in arcsec
       [(-500.0, 100.0), (400.0, 200.0)]>
   >>> c[0]
   <SkyCoord (HelioProjective: D0=149597870.7 km, obstime=None, L0=0.0 deg, B0=0.0 deg, rsun=695508.0 km): (Tx, Ty) in arcsec
       (-500.0, 100.0)>


Accessing Coordinates
^^^^^^^^^^^^^^^^^^^^^

Individual coordinates can be accessed via attributes on the SkyCoord object,
but the names of the components of the coordinates for each frame differ. For a full
description of all the properties of the frames see `sunpy.coordinates.frames`.

``HelioProjective``
###################

For the helioprojective frame the coordinates are access as ``Tx`` and ``Ty`` representing theta x and y. These are the same coordinates that are often referred to as 'solar-x' and 'solar-y'.

  >>> c = SkyCoord(-500*u.arcsec, 100*u.arcsec, frame=frames.Helioprojective)
  >>> c.Tx
  <Longitude180 -500.0 arcsec>
  >>> c.Ty
  <Latitude 100.0 arcsec>

``Heliocentric``
################

Heliocentric normally a Cartesian frame so the coordinates are accessed as ``x,y,z``:

  >>> c = SkyCoord(-72241.0*u.km, 361206.1*u.km, 589951.4*u.km, frame=frames.Heliocentric)
  >>> c.x
  <Quantity -72241.0 km>
  >>> c.y
  <Quantity 361206.1 km>
  >>> c.z
  <Quantity 589951.4 km>

``HeliographicStonyhurst`` and ``HeliographicCarrington``
#########################################################

Both the heliographic frames use latitude, longitude and radius which are accessed as follows:

   >>> c = SkyCoord(70*u.deg, -30*u.deg, frame=frames.HeliographicStonyhurst)
   >>> c.lat
   <Latitude -30.0 deg>
   >>> c.lon
   <Longitude180 70.0 deg>
   >>> c.radius
   <Distance 695508.0 km>


Observer Location Information
-----------------------------

Both ``Helioprojective`` and ``Heliocentric`` frames are defined by the location of the observer. For example in ``Helioprojective`` the observer is at the origin of the coordinate system. This information is encoded in the ``Helioprojecitve`` and ``Heliocentric`` frames as the attributes `L0`, `B0` and `D0` which are heliographic longitude, latitude and distance from the center of the Sun. These attributes are automatically populated from meta data when coordinate frames are created using map.

It is possible to convert from a ``Helioprojective`` frame with one observer location to another ``Helioprojecitve`` frame with a different observer location, by converting through ``Heliographic``, this does involve making an assumption of the radius of the Sun to calculate the position on the solar sphere. The conversion can be performed as follows::

  # Input coordinate
  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective)
  # Define the location of the new observer as a Helioprojective frame
  >>> hpc_out = sunpy.coordinates.Helioprojective(L0=10*u.deg)
  # Perform the conversion
  >>> hpc2 = hpc1.transform_to(hpc_out)

An example with two maps, i.e. ``aia`` and ``stereo``::

  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=aia.coordinate_frame)
  >>> hpc2 = hpc1.transform_to(stereo.coordinate_frame)


Design of the Coordinates Module
--------------------------------

This module works by defining a collection of ``Frames``
(`sunpy.coordinates.frames`), which exists on a transformation graph, where the
transformations between the coordinate frames are then defined and registered
with the transformation graph (`sunpy.coordinates.transformations`). Currently,
the SunPy frames are not transformable to the frames in Astropy, as there is no
transformation defined between the two sets of frames.

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
with a frame, some examples in `sunpy.coordinates` are ``obstime`` or ``L0`` and
``B0`` for observer location. Only the frames where this data is meaningful have
these attributes, i.e. only the Helioprojective frames have ``L0`` and ``B0``.
However, when you transform into another frame and then back to a projective
frame using `SkyCoord` it will remember the attributes previously provided, and
repopulate the final frame with them. If you were to do transformations using
the Frames alone this would not happen.

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
``heliographic_longitude``, ``heliographic_latitude`` and ``dsun``.

If you want to obtain a un-realized coordinate frame corresponding to a
`~sunpy.map.GenericMap` object you can do the following::

  >>> from astropy.wcs.utils import wcs_to_celestial_frame
  >>> import sunpy.map
  >>> from sunpy.data.sample import AIA_171_IMAGE

  >>> amap = sunpy.map.Map(AIA_171_IMAGE)

  >>> wcs_to_celestial_frame(amap.wcs)
  <Helioprojective Frame (D0=1.48940609627e+11 m, obstime=2011-03-1910:54:00.340000, L0=0.0 deg, B0=-7.064078 deg, rsun=695508.0 km)>



sunpy.coordinates Package
-------------------------

.. automodapi:: sunpy.coordinates.frames
    :headings: ^#

.. automodapi:: sunpy.coordinates.transformations
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

This package was developed by Pritish Chakraborty as part of GSOC 2014 and Stuart Mumford.
