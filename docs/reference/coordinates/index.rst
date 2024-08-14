.. _sunpy-coordinates-reference:

Coordinates (`sunpy.coordinates`)
*********************************

This sub-package contains:

* A robust framework for working with solar-physics coordinate systems
* Functions to obtain the locations of solar-system bodies (`sunpy.coordinates.ephemeris`)
* Functions to calculate Sun-specific coordinate information (`sunpy.coordinates.sun`)
* Bridge module to enable the use of the `~astropy.coordinates.SkyCoord` API to perform computations using `SPICE <https://naif.jpl.nasa.gov/naif/>`__ kernels (`sunpy.coordinates.spice`)

The SunPy coordinate framework extends the
:ref:`Astropy coordinates framework <astropy:astropy-coordinates>`.
The :ref:`coordinates guide <sunpy-topic-guide-coordinates-index>` provides in depth discussion of the structure and concepts underlying the coordinates framework.

.. _sunpy-coordinate-systems:

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
   * - Geomagnetic
     - MAG
     - `~sunpy.coordinates.frames.Geomagnetic`
     -
   * - Solar Magnetic
     - SM
     - `~sunpy.coordinates.frames.SolarMagnetic`
     -
   * - GeocentricSolarMagnetospheric
     - GSM
     - `~sunpy.coordinates.frames.GeocentricSolarMagnetospheric`
     -


For a description of these coordinate systems,
see `Thompson (2006) <https://doi.org/10.1051/0004-6361:20054262>`_
and `Franz & Harper (2002) <https://doi.org/10.1016/S0032-0633(01)00119-2>`_
(and `corrected version <https://www2.mps.mpg.de/homes/fraenz/systems/systems3art/systems3art.html>`_).

Reference/API
=============

.. automodapi:: sunpy.coordinates

.. automodapi:: sunpy.coordinates.ephemeris

.. automodapi:: sunpy.coordinates.spice

.. automodapi:: sunpy.coordinates.sun

.. automodapi:: sunpy.coordinates.utils
    :no-inheritance-diagram:

.. _sunpy-coordinates-other-api:

Reference/API for supporting coordinates modules
================================================

The parts of the following modules that are useful to a typical user are already imported into the `sunpy.coordinates` namespace.

.. automodapi:: sunpy.coordinates.frames

.. automodapi:: sunpy.coordinates.screens

.. automodapi:: sunpy.coordinates.metaframes

.. automodapi:: sunpy.coordinates.wcs_utils

Attribution
===========

Some of this documentation was adapted from Astropy under the terms of the `BSD
License
<https://raw.githubusercontent.com/astropy/astropy/master/LICENSE.rst>`_.

This sub-package was initially developed by Pritish Chakraborty as part of GSOC 2014 and Stuart Mumford.
