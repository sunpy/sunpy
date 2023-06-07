.. _sunpy-topic-guide-coordinates-index:

***********
Coordinates
***********

The `sunpy.coordinates` sub-package contains:

* A robust framework for working with solar-physics coordinate systems
* Functions to obtain the locations of solar-system bodies (`sunpy.coordinates.ephemeris`)
* Functions to calculate Sun-specific coordinate information (`sunpy.coordinates.sun`)

The sunpy coordinate framework extends the :ref:`Astropy coordinates framework <astropy:astropy-coordinates>`.

This page contains an overview of how coordinates work in sunpy.
See the following pages for detailed discussion of specific topics:

.. toctree::
   :maxdepth: 1

   carrington
   rotatedsunframe
   velocities
   wcs
   rsun

Design of the Coordinates Sub-Package
=====================================

This sub-package works by defining a collection of "Frames" (`sunpy.coordinates.frames`), which exists on a transformation graph, where the transformations between the coordinate frames are then defined and registered with the transformation graph (`sunpy.coordinates`).
It is also possible to transform SunPy frames to Astropy frames.

Positions within these "Frames" are stored as a "Representation" of a coordinate, a representation being a description of a point in a Cartesian, spherical or cylindrical system (see :ref:`astropy-coordinates-representations`).
A frame that contains a representation of one or many points is said to have been 'realized'.

For a more in depth look at the design and concepts of the Astropy coordinates system see :ref:`astropy-coordinates-overview`

Frames and SkyCoord
-------------------

The `~astropy.coordinates.SkyCoord` class is a high level wrapper around the `astropy.coordinates` sub-package.
It provides an easier way to create and transform coordinates, by using string representations for frames rather than the classes themselves and some other usability improvements.
For more information see the `~astropy.coordinates.SkyCoord` documentation.

The main advantage provided by `~astropy.coordinates.SkyCoord` is the support it provides for caching Frame attributes.
Frame attributes are extra data specified with a frame, some examples in `sunpy.coordinates` are ``obstime`` or ``observer`` for observer location.
Only the frames where this data is meaningful have these attributes, i.e., only the Helioprojective frames have ``observer``.
However, when you transform into another frame and then back to a projective frame using `~astropy.coordinates.SkyCoord` it will remember the attributes previously provided, and repopulate the final frame with them.
If you were to do transformations using the Frames alone this would not happen.
