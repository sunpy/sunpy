"""
The WCS package provides functions to parse World Coordinate System (WCS)
coordinates for solar images as well as convert between various solar
coordinate systems. The solar coordinates supported are

* Helioprojective-Cartesian (HPC): The most often used solar coordinate
    system. Describes positions on the Sun as angles measured from the
    center of the solar disk (usually in arcseconds) using cartesian
    coordinates (X, Y)
* Helioprojective-Radial (HPR): Describes positions on the Sun using angles,
    similar to HPC, but uses a radial coordinate (rho, psi) system centered
    on solar disk where psi is measured in the counter clock wise direction.
* Heliocentric-Cartesian (HCC): The same as HPC but with positions expressed
    in true (deprojected) physical distances instead of angles on the
    celestial sphere.
* Heliocentric-Radial (HCR): The same as HPR but with rho expressed in
    true (deprojected) physical distances instead of angles on the celestial
    sphere.
* Stonyhurst-Heliographic (HG): Expressed positions on the Sun using
    longitude and latitude on the solar sphere but with the origin which is
    at the intersection of the solar equator and the central meridian as
    seen from Earth. This means that the coordinate system remains fixed
    with respect to Earth while the Sun rotates underneath it.
* Carrington-Heliographic (HG): Carrington longitude is offset
    from Stonyhurst longitude by a time-dependent scalar value, L0. At the
    start of each Carrington rotation, L0 = 360, and steadily decreases
    until it reaches L0 = 0, at which point the next Carrington rotation
    starts.

Some definitions

* b0: Tilt of the solar North rotational axis toward the observer
    (helio- graphic latitude of the observer). Note that SOLAR_B0,
    HGLT_OBS, and CRLT_OBS are all synonyms.
* l0: Carrington longitude of central meridian as seen from Earth.
* dsun_meters: Distance between observer and the Sun. Default is 1 AU.
* rsun_meters: Radius of the Sun in meters. Default is 6.955e8 meters. This valued is stored
  locally in this module and can be modified if necessary.

References
----------
| Thompson (2006), A&A, 449, 791 <http://dx.doi.org/10.1051/0004-6361:20054262>
| PDF <http://fits.gsfc.nasa.gov/wcs/coordinates.pdf>
"""
from __future__ import absolute_import

from sunpy.wcs.wcs import *
