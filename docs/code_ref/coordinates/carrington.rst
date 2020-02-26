.. _sunpy-coordinates-carrington:

Calculating Carrington longitude
================================

Carrington coordinates are a heliographic coordinate system, where the longitude is determined assuming that the Sun has rotated at a sidereal period of ~25.38 days starting at zero longitude at a specific time on 1853 November 9.
For a given observer at a given observation time, the center of the disk of the Sun as seen by that observer (also known as the apparent sub-observer point) has a Carrington latitude and longitude.
When the observer is at Earth, the apparent sub-Earth Carrington latitude and longitude are also known as |B0| and |L0|, respectively.

Quick summary
-------------
SunPy calculates the sub-observer Carrington longitude in a manner that is consistent with |AA| for Earth, but enables co-alignment of images from different observatories at different locations/velocities in the solar system.

Observer effects
----------------
There are two primary effects for why the apparent sub-observer point is not the same as the "true" sub-observer point.

Light travel time
^^^^^^^^^^^^^^^^^

It takes time for light to travel from the Sun to any observer (e.g., ~500 seconds for an observer at Earth).
Thus, the farther the observer is from the Sun, the farther back in time the Sun will appear to be.
This effect is sometimes called "planetary aberration".

As an additional detail, the observer is closer to the sub-observer point than to the center of the Sun (i.e., the radius of the Sun), which corresponds to ~2.3 lightseconds.

Observer motion
^^^^^^^^^^^^^^^

The motion of an observer will shift the apparent location of objects, and this effect is called "stellar aberration".
A more subtle question is whether observer motion also shifts the apparent sub-observer Carrington longitude.

The methodology of the |AA| is that stellar aberration causes the apparent ecliptic longitude of an observer to be different from its true ecliptic longitude.
The methodology then uses this apparent ecliptic longitude to calculate the sub-observer Carrington longitude.
This results in a shift of ~0.006 degrees (~20 arcseconds) in Carrington longitude.

However, shifting the sub-observer Carrington longitude due to observer motion is not always appropriate.
For certain types of observations (e.g., imagery of features on the surface of the Sun), stellar aberration does not shift the positions of solar features **relative** to the center of the disk of the Sun; everything shifts in tandem.
Thus, for the purposes of co-aligning images from different observatories, it is critical to **exclude** any correction for observer motion.

Background information
----------------------

The IAU definition
^^^^^^^^^^^^^^^^^^

Since the original definition of Carrington coordinates [#Carrington]_, there have been a number of refinements of the definition, but there have been inconsistencies with how those definitions are interpreted.
Possible points of contention include the reference epoch, the rotation period, and the handling of observer effects.

The most recent definition was essentially laid down in 2007 by the IAU [#IAU]_:

* The reference epoch is J2000.0 TDB (2000 January 1 12:00 TDB).
  The difference between Barycentric Dynamical Time (TDB) and Terrestrial Time (TT) is very small, but both are appreciably different from Coordinated Universal Time (UTC).
* The longitude of the prime meridian (|W0|) at the reference epoch is 84.176 degrees.
  This longitude is the "true" longitude, without accounting for any effects that would modify the apparent prime meridian for an observer (e.g., light travel time).
* The sidereal rotation rate is exactly 14.1844 degrees per day.
  This is close to, but not exactly the same as, a sidereal rotation period of 25.38 days or a mean synodic period of 27.2753 days.
* The ICRS coordinates of the Sun's north pole of rotation is given as a right ascension of 286.13 degrees and a declination of 63.87 degrees.

Using the IAU definition correctly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Importantly, the IAU parameters (specifically, |W0|) were tuned following an investigation by |AA| [#AA]_.
The prescription of |AA| included both a correction for light travel time and a correction for stellar aberration due to Earth's motion.
When using this prescription, |AA| found that a |W0| value of 84.176 degrees minimized the differences between the modern approach and the earlier approach.
Thus, it is critical to follow the prescription of |AA| when using the IAU parameters, otherwise the results will not match the published values in |AA|.

The approach in SunPy
---------------------

SunPy determines apparent sub-observer Carrington longitude as a **relative** calculation to the sub-Earth Carrington longitude (|L0|).
This approach includes a correction for the difference in light travel time, but **excludes** any correction for stellar aberration due to observer motion.
The exclusion of a stellar-aberration correction is appropriate for purposes such as image co-alignment, where Carrington longitudes must be consistently associated with features on the Sun's surface.

However, calculating |L0| from the current IAU parameters requires following the manner prescribed by |AA|.
This manner includes both a correction for light travel time and a correction for stellar aberration due to Earth's motion.
By calculating |L0| this way, SunPy matches the values published by |AA| and minimizes errors with historical definitions.

Comparisons to other sources
----------------------------

Compared to |AA|
^^^^^^^^^^^^^^^^

|AA| publishes the apparent sub-Earth Carrington longitude (|L0|) to a precision of 0.01 degrees.
|AA| does not appear to account for the difference in light travel time between the sub-Earth point and the center of the Sun (~2.3 lightseconds), which results in a fixed discrepancy of ~1.4 arcseconds in |L0|.
This small discrepancy occasionally results in disagreement at the printed precision due to rounding.

Compared to SPICE
^^^^^^^^^^^^^^^^^

In `SPICE <https://naif.jpl.nasa.gov/naif/>`_, the apparent sub-observer Carrington longitude can be calculated using the function `subpnt <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html>`_.
`subpnt`_ allows for two types of aberration correction through the input argument ``abcorr``:

* "LT", which is just the correction for light travel time
* "LT+S", which is the correction for light travel time ("LT") plus the correction for observer motion ("S")

SunPy calculates the apparent sub-Earth Carrington longitude (|L0|) in a manner equivalent to specifying "LT+S" for Earth.
Recall that SunPy includes the correction for observer motion to correctly use the IAU definition.
The discrepancy is no greater than 0.2 arcseconds.

SunPy calculates the apparent sub-observer Carrington longitude in a manner equivalent to specifying "LT" for the observer and then adding the difference between "LT" and "LT+S" for Earth.
The discrepancy is no greater than 0.2 arcseconds.

Compared to JPL Horizons
^^^^^^^^^^^^^^^^^^^^^^^^

In `JPL Horizons <https://ssd.jpl.nasa.gov/?horizons>`_, one can request the "Obs sub-long & sub-lat".
JPL Horizons appears to start from the IAU parameters and to include the correction for light travel time but not the correction for observer motion (i.e., equivalent to specifying "LT" to `subpnt`_ in `SPICE`_).
Thus, JPL Horizons reports values for sub-Earth Carrington longitude (|L0|) that are ~20 arcseconds offset from SunPy (as well as |AA|)

Compared to SunSPICE
^^^^^^^^^^^^^^^^^^^^

In `SunSPICE <https://stereo-ssc.nascom.nasa.gov/sunspice.shtml>`_, one can convert to and from the "Carrington" coordinate system using the function `convert_sunspice_coord <https://hesperia.gsfc.nasa.gov/ssw/packages/sunspice/idl/convert_sunspice_coord.pro>`_.
However, these Carrington longitudes are "true" rather than "apparent" because the observer is not specified, so there are no corrections for light travel time or for observer motion.
Thus, the discrepancy is ~0.82 degrees.

Footnotes
---------

.. [#Carrington] Carrington (1863), *Observations of the Spots on the Sun*, p. 244
.. [#IAU] Seidelmann et al. (2007), "Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2006", `<http://dx.doi.org/10.1007/s10569-007-9072-y>`__
.. [#AA] Urban & Kaplan (2007), "Investigation of Change in the Computational Technique of the Sun’s Physical Ephemeris in The Astronomical Almanac", `<http://asa.hmnao.com/static/files/sun_rotation_change.pdf>`__

.. |AA| replace:: *The Astronomical Almanac*
.. |B0| replace:: B\ :sub:`0`
.. |L0| replace:: L\ :sub:`0`
.. |W0| replace:: W\ :sub:`0`
