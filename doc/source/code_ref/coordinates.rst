SunPy coordinates
=================

The SunPy coordinates submodule is an implementation of the common solar physics
coordinate frames in to the :ref:`Astropy coordinates system <astropy:astropy-coordinates>`.

.. warning::

   The accuracy of the transformations in this module have not been rigorously
   verified. They have been compared to `sunpy.wcs` and match to numerical
   precision. Independent verification will be added at a later date.


Getting Started
---------------

The easiest interface to the coordinates module is through the `~astropy.coordinates.SkyCoord` class::

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord
  >>> import sunpy.coordinates
  >>> c = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame='helioprojective')
  >>> c = SkyCoord(90*u.deg, -76*u.deg, frame='helioprojectiveradial')
  >>> c = SkyCoord(x=-72241.0*u.km, y=361206.1*u.km, z=589951.4*u.km, frame='heliocentric')
  >>> c = SkyCoord(70*u.deg, -30*u.deg, frame='heliographicstonyhurst')
  >>> c
  <SkyCoord (HelioGraphicStonyhurst: dateobs=None): (lon, lat, rad) in (deg, deg, km)
      (70.0, -30.0, 695508.0)>


SunPy implements supports for the following solar physics coordinate systems:

* Helioprojective (Cartesian) `~sunpy.coordinates.frames.HelioProjective`
* Helioprojective (Radial) `~sunpy.coordinates.frames.HelioProjectiveRadial`
* Heliocentric `~sunpy.coordinates.frames.HelioCentric`
* Heliographic Stonyhurst `~sunpy.coordinates.frames.HelioGraphicStonyhurst`
* Heliographic Carrington `~sunpy.coordinates.frames.HelioGraphicCarrington`

for a complete description of these frames, see `sunpy.coordinates.frames`.


`~astropy.coordinates.SkyCoord` and all other `~astropy.coordinates` objects
also support array coordinates. These work the same as single-value coordinates,
but they store multiple coordinates in a single object. When you're going to
apply the same operation to many different coordinates, this is a better choice
than a list of `~astropy.coordinates.SkyCoord` objects, because it will be
*much* faster than applying the operation to each
`~astropy.coordinates.SkyCoord` in a for loop.
::

   >>> c = SkyCoord([-500, 400]*u.arcsec, [100, 200]*u.arcsec, frame='helioprojective')
   >>> c
   <SkyCoord (HelioProjective: D0=149597870.7 km, dateobs=None, L0=0.0 deg, B0=0.0 deg, RSun=695508.0 km): (Tx, Ty) in arcsec
       [(-500.0, 100.0), (400.0, 200.0)]>
   >>> c[0]
   <SkyCoord (HelioProjective: D0=149597870.7 km, dateobs=None, L0=0.0 deg, B0=0.0 deg, RSun=695508.0 km): (Tx, Ty) in arcsec
       (-500.0, 100.0)>


Accessing Coordinates
^^^^^^^^^^^^^^^^^^^^^

Individual coordinates can be accessed via attributes on the SkyCoord object,
the names of the components of the coordinates for each frame differ. For a full
description of all the properties of the frames see `sunpy.coordinates.frames`.

``HelioProjective``::

  >>> c = SkyCoord(-500*u.arcsec, 100*u.arcsec, frame='helioprojective')
  >>> c.Tx
  <Longitude180 -500.0 arcsec>
  >>> c.Ty
  <Latitude 100.0 arcsec>

``Heliocentric``::

  >>> c = SkyCoord(-72241.0*u.km, 361206.1*u.km, 589951.4*u.km, frame='heliocentric')
  >>> c.x
  <Quantity -72241.0 km>
  >>> c.y
  <Quantity 361206.1 km>
  >>> c.z
  <Quantity 589951.4 km>

``HeliographicStonyhurst``::

   >>> c = SkyCoord(70*u.deg, -30*u.deg, frame='heliographicstonyhurst')
   >>> c.lat
   <Latitude -30.0 deg>
   >>> c.lon
   <Longitude180 70.0 deg>
   >>> c.rad
   <Distance 695508.0 km>


Design of the Coordinates Module
--------------------------------

This module works by defining a collection of ``Frames``
(`sunpy.coordinates.frames`), which exists on a transformation graph, where the
transformations between the coordinate frames are then defined and registered
with the transformation graph (`sunpy.coordinates.transformations`).

Positions within these ``Frames`` are stored as a ``Representation`` of a
coordinate, a representation is a description of a point in a Cartesian,
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
with a frame, some examples in `sunpy.coordinates` are ``dateobs`` or ``L0`` and
``B0`` for observer location. Only the frames where this data is meaningful have
these attributes, i.e. only the Helioprojective frames have ``L0`` and ``B0``.
However, when you transform into another frame and then back to a projective
frame using `SkyCoord` it will remember the attributes previously provided, and
repopulate the final frame with them. If you were to do transformations using
the Frames alone this would not happen.

The most important implication for this in `sunpy.coordinates` is the ``RSun``
parameter in the projective frames. If you create a projective frame with a
``RSun`` attribute, if you convert back to a projective frame it will be set
correctly. It should also be noted that, if you create a Heliographic frame and
then transform to a projective frame with an ``RSun`` attribute, it will not
match the ``rad`` coordinate in the Heliographic frame. This is because you may
mean to be describing a point above the defined 'surface' of the Sun.


Coordinates and WCS
-------------------

The `sunpy.coordinates` package provides a mapping between FITS-WCS CTYPE
convention and the coordinate frames as defined in `sunpy.coordinates`. This is
provided in the `sunpy.coordinates.wcs_utils` module and adds the mappings to
the `astropy.wcs.utils.WCS_FRAME_MAPPINGS` list. This list is used by packages
such as ``wcsaxes`` to convert from `astropy.wcs.WCS` objects to coordinate
frames.

The `sunpy.map.GenricMap` class creates `astropy.wcs.WCS` objects as
``amap.wcs``, however, it adds some extra attributes to the `~astropy.wcs.WCS`
object to be able to fully specify the coordinate frame. It adds
``heliographic_longitude``, ``heliographic_latitude`` and ``dsun``.

If you want to obtain a un-realized coordinate frame corresponding to a
`~sunpy.map.GenericMap` object you can do the following::

  >>> from astropy.wcs.utils import wcs_to_celestial_frame
  >>> import sunpy.coordinates
  >>> import sunpy.map
  >>> from sunpy.data.sample import AIA_171_IMAGE

  >>> amap = sunpy.map.Map(AIA_171_IMAGE)

  >>> wcs_to_celestial_frame(amap.wcs)
  <Helioprojective Frame (D0=1.48940609627e+11 m, dateobs=2011-03-1910:54:00.340000, L0=0.0 deg, B0=-7.064078 deg, RSun=695508.0 km)>



sunpy.coordinates Package
-------------------------

.. automodapi:: sunpy.coordinates.frames
    :headings: ^#

.. automodapi:: sunpy.coordinates.transformations
    :headings: ^#

.. automodapi:: sunpy.coordinates.representation
    :headings: ^#

.. automodapi:: sunpy.coordinates.wcs_utils
    :headings: ^#


Attribution
-----------

Some of this documentation was borrowed from Astropy under the terms of the `BSD
License
<https://raw.githubusercontent.com/astropy/astropy/master/licenses/LICENSE.rst>`_.

This package was developed by Pritish Chakraborty as part of GSOC 2014.
