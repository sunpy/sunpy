SunPy coordinates
=================

The SunPy coordinates submodule is an implementation of the common solar physics
coordinate frames in to the :ref:`Astropy coordinates system <astropy:astropy-coordinates>`.

.. warning::

    This submodule is in beta, the accuracy of the transformations have not 
    been verified. See `sunpy.wcs` for the previous version of these 
    transformations.


Getting Started
---------------

The easiest interface to the coordinates module is through the `~astropy.coordinates.SkyCoord` module::

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord
  >>> import sunpy.coordinates
  >>> c = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame='helioprojective')
  >>> c = SkyCoord(90*u.deg, -76*u.deg, frame='helioprojectiveradial')
  >>> c = SkyCoord(x=-72241.0*u.km, y=361206.1*u.km, z=589951.4*u.km, frame='heliocentric')
  >>> c = SkyCoord(70*u.deg, -30*u.deg, frame='heliographicstonyhurst')
  >>> c
  <SkyCoord (HelioGraphicStonyhurst: dateobs=None, RSun=695508.0 km): (lon, lat, rad) in (deg, deg, km)
      (70.0, -30.0, 695508.0)>

`~astropy.coordinates.SkyCoord` and all other `~astropy.coordinates` objects
also support array coordinates. These work the same as single-value coordinates,
but they store multiple coordinates in a single object. When you're going to
apply the same operation to many different coordinates (say, from a catalog),
this is a better choice than a list of `~astropy.coordinates.SkyCoord` objects,
because it will be *much* faster than applying the operation to each
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
<https://raw.githubusercontent.com/astropy/astropy/master/licenses/LICENSE.rst>`_
