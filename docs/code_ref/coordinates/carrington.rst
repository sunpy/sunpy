.. _sunpy-coordinates-carrington:

Calculating Carrington longitude
================================

Carrington coordinates are a heliographic coordinate system, where the longitude is determined assuming that the Sun has rotated at a sidereal period of ~25.38 days starting at zero longitude at a specific time on 1853 November 9.
When working in Carrington coordinates, the primary coordinate of reference is the Carrington latitude and longitude of the center of the disk of the Sun as seen by an observer at Earth at a given observation time.
This point is the apparent sub-Earth point, and this latitude and longitude are also known as |B0| and |L0|, respectively.

Quick summary
-------------

SunPy calculates the apparent sub-Earth Carrington longitude (|L0|) from the current IAU parameters in the manner prescribed by |AA|.
This manner includes both a correction for light travel time and a correction for aberration due to Earth's motion in order to match the values published by |AA| and to minimize errors with historical definitions.

For a given observer, SunPy calculates the apparent sub-observer Carrington longitude as a **relative** calculation to |L0|.
This approach includes a correction for the difference in light travel time, but **excludes** any correction for aberration due to observer motion, in order to enable the consistent association of Carrington longitudes to solar features.

The IAU definition
------------------

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

Observer effects
----------------

Light travel time
^^^^^^^^^^^^^^^^^

It takes time for light to travel from the Sun to any observer (e.g., ~500 seconds for an observer at Earth).
Thus, the farther the observer is from the Sun, the farther back in time the Sun will appear to be.

Observer motion
^^^^^^^^^^^^^^^

The motion of an observer will shift the apparent location of objects, and this effect is called "aberration".
A more subtle question is whether observer motion also shifts the apparent sub-observer Carrington longitude.

The methodology of the |AA| is that aberration causes the apparent ecliptic longitude of an observer to be different from its true ecliptic longitude.
The methodology then uses this apparent ecliptic longitude to calculate the sub-observer Carrington longitude.
This results in a shift of ~0.006 degrees in Carrington longitude.

However, shifting the sub-observer Carrington longitude due to observer motion may not be appropriate.

Using the IAU definition correctly
----------------------------------

Importantly, the IAU parameters (specifically, |W0|) were tuned following an investigation by |AA| [#AA]_.
The prescription of |AA| included both a correction for light travel time and a correction for aberration due to Earth's motion [#CSPICE1]_.
When using this prescription, |AA| found that a |W0| value of 84.176 degrees minimized the differences between the modern approach and the earlier approach.
Thus, it is critical to follow the prescription of |AA| when using the IAU parameters, otherwise the results will not match the published values in |AA|.

If one applies the correction for light travel time, but *not* the correction for aberration due to Earth's motion [#CSPICE2]_, the calculated apparent sub-observer Carrington longitudes will be shifted by ~0.006 degrees.
This is the methodology of JPL Horizons [#Horizons]_, and that is the reason for discrepancies between their numbers and |AA|.

Ignoring observer motion
------------------------

SunPy does not shift the sub-observer Carrington longitude due to observer motion.
Even though aberration shifts the apparent angular direction of the Sun, aberration does not shift the positions of solar features **relative** to the center of the disk of the Sun.
For example, for a given observer location, a sunspot that appears to be at the center of the disk will always be at the center of the disk irrespective of the amount of observer motion.
Thus, for the purposes of co-aligning images from different observatories, it is critical to **exclude** any correction for aberration due to observer motion.


Additional considerations
-------------------------
The SunPy methodology does not, by default, account for the fact that the sub-observer point is about ~2 light-seconds closer to the observer than the center of the Sun.
We believe that this correction is not part of the |AA| prescription.
This correction corresponds to ~1.4 arcseconds in Carrington longitude.
The function :func:`sunpy.coordinates.sun.L0` has the optional keyword ``nearest_point`` to include this correction.
When including this correction, SunPy currently agrees with CSPICE calculations of apparent Carrington longitude to ~0.2 arcseconds.

.. |AA| replace:: *The Astronomical Almanac*
.. |B0| replace:: B\ :sub:`0`
.. |L0| replace:: L\ :sub:`0`
.. |W0| replace:: W\ :sub:`0`

.. [#Carrington] Carrington (1863), *Observations of the Spots on the Sun*, p. 244
.. [#IAU] Seidelmann et al. (2007), "Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2006", `<http://dx.doi.org/10.1007/s10569-007-9072-y>`__
.. [#AA] Urban & Kaplan (2007), "Investigation of Change in the Computational Technique of the Sun’s Physical Ephemeris in The Astronomical Almanac", `<http://asa.hmnao.com/static/files/sun_rotation_change.pdf>`__
.. [#CSPICE1] In CSPICE, this prescription is equivalent to specifying "LT+S" for ``abcorr`` to the function `subpnt <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html>`_.
.. [#CSPICE2] In CSPICE, this difference is equivalent to specifying "LT" instead of "LT+S" for ``abcorr`` to the function `subpnt`_.
.. [#Horizons] `<https://ssd.jpl.nasa.gov/?horizons>`__
