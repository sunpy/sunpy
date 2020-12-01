"""
Sun-specific coordinate calculations
"""
import numpy as np

import astropy.units as u
from astropy.constants import c as speed_of_light
from astropy.coordinates import (
    ITRS,
    AltAz,
    Angle,
    Distance,
    GeocentricMeanEcliptic,
    GeocentricTrueEcliptic,
    HeliocentricMeanEcliptic,
    Latitude,
    Longitude,
    SkyCoord,
    get_body_barycentric,
)
from astropy.coordinates.builtin_frames.utils import get_jd12
from astropy.coordinates.representation import CartesianRepresentation, SphericalRepresentation
# Import erfa via astropy to make sure we are using the same ERFA library as Astropy
from astropy.coordinates.sky_coordinate import erfa
from astropy.time import Time

from sunpy import log
from sunpy.sun import constants
from sunpy.time import parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring
from .ephemeris import get_earth
from .frames import HeliographicStonyhurst
from .transformations import _SOLAR_NORTH_POLE_HCRS, _SUN_DETILT_MATRIX

__author__ = "Albert Y. Shih"
__email__ = "ayshih@gmail.com"

__all__ = [
    "angular_radius", "sky_position", "carrington_rotation_number",
    "carrington_rotation_time",
    "true_longitude", "apparent_longitude", "true_latitude", "apparent_latitude",
    "mean_obliquity_of_ecliptic", "true_rightascension", "true_declination",
    "true_obliquity_of_ecliptic", "apparent_rightascension", "apparent_declination",
    "print_params",
    "B0", "L0", "P", "earth_distance", "orientation"
]


@add_common_docstring(**_variables_for_parse_time_docstring())
def angular_radius(t='now'):
    """
    Return the angular radius of the Sun as viewed from Earth.

    The tangent vector from the Earth to the edge of the Sun forms a
    right-angle triangle with the radius of the Sun as the far side and the
    Sun-Earth distance as the hypotenuse.  Thus, the sine of the angular
    radius of the Sun is ratio of these two distances.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format
    """
    return _angular_radius(constants.radius, earth_distance(t))


def _angular_radius(sol_radius, distance):
    solar_semidiameter_rad = np.arcsin(sol_radius / distance)
    return Angle(solar_semidiameter_rad.to(u.arcsec))


@add_common_docstring(**_variables_for_parse_time_docstring())
def sky_position(t='now', equinox_of_date=True):
    """
    Returns the apparent position of the Sun (right ascension and declination) on the
    celestial sphere using the equatorial coordinate system, referred to the true equinox of date
    (as default).  Corrections for nutation and aberration (for Earth motion) are included.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format

    equinox_of_date : `bool`
        If True, output is referred to the true equinox of date.  Otherwise, output is referred to
        the J2000.0 epoch (ICRF orientation, not dynamical orientation).
    """
    ra = apparent_rightascension(t, equinox_of_date=equinox_of_date)
    dec = apparent_declination(t, equinox_of_date=equinox_of_date)
    return ra, dec


def carrington_rotation_time(crot):
    """
    Return the time of a given Carrington rotation.

    The round-trip from this method to `carrington_rotation_number` has
    absolute errors of < 0.11 seconds.

    Parameters
    ----------
    crot : float
        Carrington rotation number. Can be a fractional cycle number.

    Returns
    -------
    astropy.time.Time
    """
    estimate = (constants.mean_synodic_period.to_value('day') *
                (crot - 1)) + constants.first_carrington_rotation

    # The above estimate is inaccurate (see comments below in carrington_rotation_number),
    # so put the estimate into carrington_rotation_number to determine a correction amount
    def refine(estimate):
        crot_estimate = carrington_rotation_number(estimate)
        dcrot = crot - crot_estimate
        # Correct the estimate using a linear fraction of the Carrington rotation period
        return estimate + (dcrot * constants.mean_synodic_period)

    # Perform two iterations of the correction to achieve sub-second accuracy
    estimate = refine(estimate)
    estimate = refine(estimate)
    return Time(estimate, scale='tt', format='jd')


@add_common_docstring(**_variables_for_parse_time_docstring())
def carrington_rotation_number(t='now'):
    """
    Return the Carrington rotation number.  Each whole rotation number marks when the Sun's prime
    meridian coincides with the central meridian as seen from Earth, with the first rotation
    starting on 1853 November 9.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format
    """
    time = parse_time(t)

    # Estimate the Carrington rotation number by dividing the time that has elapsed since
    # JD 2398167.4 (late in the day on 1853 Nov 9), see Astronomical Algorithms (Meeus 1998, p.191),
    # by the mean synodic period (27.2753 days)
    estimate = (time - constants.first_carrington_rotation) / constants.mean_synodic_period + 1
    estimate_int, estimate_frac = divmod(estimate, 1)

    # The fractional rotation number from the above estimate is inaccurate, so calculate the actual
    # fractional rotation number from the longitude of the central meridian (L0)
    actual_frac = 1 - L0(time).to('deg').value / 360

    # Calculate any adjustment to the integer rotation number due to wrapping
    wrap_adjustment = np.around(estimate_frac - actual_frac)

    actual = estimate_int + actual_frac + wrap_adjustment

    log.debug(f"Carrington rotation number: estimate is {estimate}, actual is {actual}")

    return actual.to_value(u.one)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_longitude(t='now'):
    """
    Returns the Sun's true geometric longitude, referred to the mean equinox of date.  No
    corrections for nutation or aberration are included.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format
    """
    time = parse_time(t)

    # Calculate Earth's true geometric longitude and add 180 degrees for the Sun's longitude.
    # This approach is used because Astropy's GeocentricMeanEcliptic includes aberration.
    earth = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='gcrs', obstime=time)
    coord = earth.transform_to(HeliocentricMeanEcliptic(equinox=time))
    lon = coord.lon + 180*u.deg

    return Longitude(lon)


@add_common_docstring(**_variables_for_parse_time_docstring())
def apparent_longitude(t='now'):
    """
    Returns the Sun's apparent longitude, referred to the true equinox of date.  Corrections for
    nutation and aberration (for Earth motion) are included.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format

    Notes
    -----
    The nutation model is IAU 2000A nutation with adjustments to match IAU 2006 precession.
    """
    time = parse_time(t)
    sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='hcrs', obstime=time)
    coord = sun.transform_to(GeocentricTrueEcliptic(equinox=time))

    # Astropy's GeocentricTrueEcliptic includes both aberration and nutation
    lon = coord.lon

    return Longitude(lon)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_latitude(t='now'):
    """
    Returns the Sun's true geometric latitude, referred to the mean equinox of date.  No
    corrections for nutation or aberration are included.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format
    """
    time = parse_time(t)
    sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='hcrs', obstime=time)
    coord = sun.transform_to(GeocentricMeanEcliptic(equinox=time))

    # Astropy's GeocentricMeanEcliptic includes aberration from Earth motion, but the contribution
    # is negligible
    lat = coord.lat

    return Latitude(lat)


@add_common_docstring(**_variables_for_parse_time_docstring())
def apparent_latitude(t='now'):
    """
    Returns the Sun's apparent latitude, referred to the true equinox of date.  Corrections for
    nutation and aberration (for Earth motion) are included.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format
    """
    time = parse_time(t)
    sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='hcrs', obstime=time)
    coord = sun.transform_to(GeocentricTrueEcliptic(equinox=time))

    # Astropy's GeocentricTrueEcliptic includes both aberration and nutation
    lat = coord.lat

    return Latitude(lat)


@add_common_docstring(**_variables_for_parse_time_docstring())
def mean_obliquity_of_ecliptic(t='now'):
    """
    Returns the mean obliquity of the ecliptic, using the IAU 2006 definition.  No correction for
    nutation is included.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format
    """
    time = parse_time(t)
    jd1, jd2 = get_jd12(time, 'tt')
    obl = erfa.obl06(jd1, jd2)*u.radian
    return Angle(obl, u.arcsec)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_rightascension(t='now', equinox_of_date=True):
    """
    Returns the Sun's true geometric right ascension relative to Earth, referred to the mean equinox
    of date (as default).  No corrections for nutation or aberration are included.  The correction
    due to light travel time would be negligible, so the output is also the astrometric right
    ascension.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format

    equinox_of_date : `bool`
        If True, output is referred to the mean equinox of date.  Otherwise, output is referred to
        the J2000.0 epoch (ICRF orientation, not dynamical orientation).
    """
    if equinox_of_date:
        # Mean equinox of date
        obl = mean_obliquity_of_ecliptic(t)  # excludes nutation
        lon = true_longitude(t)
        lat = true_latitude(t)

        # See Astronomical Algorithms (Meeus 1998 p.93)
        y = np.sin(lon) * np.cos(obl) - np.tan(lat) * np.sin(obl)
        x = np.cos(lon)
        result = np.arctan2(y, x)
    else:
        # J2000.0 epoch
        # Calculate Earth's true geometric right ascension relative to the Sun and add 180 degrees.
        # This approach is used because Astropy's GCRS includes aberration.
        earth = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='gcrs', obstime=parse_time(t))
        result = earth.hcrs.ra + 180*u.deg

    return Longitude(result, u.hourangle)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_declination(t='now', equinox_of_date=True):
    """
    Returns the Sun's true geometric declination relative to Earth, referred to the mean equinox
    of date (as default).  No corrections for nutation or aberration are included.  The correction
    due to light travel time would be negligible, so the output is also the astrometric declination.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format

    equinox_of_date : `bool`
        If True, output is referred to the mean equinox of date.  Otherwise, output is referred to
        the J2000.0 epoch (ICRF orientation, not dynamical orientation).
    """
    if equinox_of_date:
        # Mean equinox of date
        obl = mean_obliquity_of_ecliptic(t)  # excludes nutation
        lon = true_longitude(t)
        lat = true_latitude(t)

        # See Astronomical Algorithms (Meeus 1998 p.93)
        result = np.arcsin(np.sin(lat) * np.cos(obl) + np.cos(lat) * np.sin(obl) * np.sin(lon))
    else:
        # J2000.0 epoch
        # Calculate Earth's true geometric declination relative to the Sun and multipy by -1.
        # This approach is used because Astropy's GCRS includes aberration.
        earth = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='gcrs', obstime=parse_time(t))
        result = -earth.hcrs.dec

    return Latitude(result, u.deg)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_obliquity_of_ecliptic(t='now'):
    """
    Returns the true obliquity of the ecliptic, using the IAU 2006 definition.  Correction for
    nutation is included.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format

    Notes
    -----
    The nutation model is IAU 2000A nutation with adjustments to match IAU 2006 precession.
    """
    time = parse_time(t)
    jd1, jd2 = get_jd12(time, 'tt')
    obl = erfa.obl06(jd1, jd2)*u.radian
    _, nut_obl = erfa.nut06a(jd1, jd2)*u.radian
    obl += nut_obl
    return Angle(obl, u.arcsec)


@add_common_docstring(**_variables_for_parse_time_docstring())
def apparent_rightascension(t='now', equinox_of_date=True):
    """
    Returns the Sun's apparent right ascension relative to Earth, referred to the true equinox
    of date (as default).  Corrections for nutation or aberration (for Earth motion) are included.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format

    equinox_of_date : `bool`
        If True, output is referred to the true equinox of date.  Otherwise, output is referred to
        the J2000.0 epoch (ICRF orientation, not dynamical orientation).
    """
    if equinox_of_date:
        # True equinox of date
        obl = true_obliquity_of_ecliptic(t)  # includes nutation
        lon = apparent_longitude(t)
        lat = apparent_latitude(t)

        # See Astronomical Algorithms (Meeus 1998 p.93)
        y = np.sin(lon) * np.cos(obl) - np.tan(lat) * np.sin(obl)
        x = np.cos(lon)
        result = np.arctan2(y, x)
    else:
        # J2000.0 epoch
        sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='hcrs', obstime=parse_time(t))
        result = sun.gcrs.ra

    return Longitude(result, u.hourangle)


@add_common_docstring(**_variables_for_parse_time_docstring())
def apparent_declination(t='now', equinox_of_date=True):
    """
    Returns the Sun's apparent declination relative to Earth, referred to the true equinox
    of date (as default).  Corrections for nutation or aberration (for Earth motion) are included.

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format

    equinox_of_date : `bool`
        If True, output is referred to the true equinox of date.  Otherwise, output is referred to
        the J2000.0 epoch (ICRF orientation, not dynamical orientation).
    """
    if equinox_of_date:
        # True equinox of date
        obl = true_obliquity_of_ecliptic(t)  # includes nutation
        lon = apparent_longitude(t)
        lat = apparent_latitude(t)

        # See Astronomical Algorithms (Meeus 1998 p.93)
        result = np.arcsin(np.sin(lat) * np.cos(obl) + np.cos(lat) * np.sin(obl) * np.sin(lon))
    else:
        # J2000.0 epoch
        sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='hcrs', obstime=parse_time(t))
        result = sun.gcrs.dec

    return Latitude(result, u.deg)


@add_common_docstring(**_variables_for_parse_time_docstring())
def print_params(t='now'):
    """
    Print out a summary of solar ephemeris. 'True' values are true geometric values referred to the
    mean equinox of date, with no corrections for nutation or aberration.  'Apparent' values are
    referred to the true equinox of date, with corrections for nutation and aberration (for Earth
    motion).

    Parameters
    ----------
    t : {parse_time_types}
        Time to use in a parse-time-compatible format
    """
    print('Solar Ephemeris for {} UTC\n'.format(parse_time(t).utc))
    print('Distance = {}'.format(earth_distance(t)))
    print('Semidiameter = {}'.format(angular_radius(t)))
    print('True (long, lat) = ({}, {})'.format(true_longitude(t).to_string(),
                                               true_latitude(t).to_string()))
    print('Apparent (long, lat) = ({}, {})'.format(apparent_longitude(t).to_string(),
                                                   apparent_latitude(t).to_string()))
    print('True (RA, Dec) = ({}, {})'.format(true_rightascension(t).to_string(),
                                             true_declination(t).to_string()))
    print('Apparent (RA, Dec) = ({}, {})'.format(apparent_rightascension(t).to_string(),
                                                 apparent_declination(t).to_string()))
    print('Heliographic long. and lat of disk center = ({}, {})'.format(L0(t).to_string(),
                                                                        B0(t).to_string()))
    print('Position angle of north pole = {}'.format(P(t)))
    print('Carrington rotation number = {}'.format(carrington_rotation_number(t)))


@add_common_docstring(**_variables_for_parse_time_docstring())
def B0(time='now'):
    """
    Return the B0 angle for the Sun at a specified time, which is the heliographic latitude of the
    of the center of the disk of the Sun as seen from Earth. The range of B0 is +/-7.23 degrees.

    Equivalent definitions include:
        * The heliographic latitude of Earth
        * The tilt of the solar North rotational axis toward Earth

    Parameters
    ----------
    time : {parse_time_types}
        Time to use in a parse_time-compatible format

    Returns
    -------
    out : `~astropy.coordinates.Latitude`
        The position angle
    """
    return Latitude(get_earth(time).lat)


# Function returns a SkyCoord's longitude in the de-tilted frame (HCRS rotated so that the Sun's
# rotation axis is aligned with the Z axis)
def _detilt_lon(coord):
    coord_detilt = coord.hcrs.cartesian.transform(_SUN_DETILT_MATRIX)
    return coord_detilt.represent_as(SphericalRepresentation).lon.to('deg')


# J2000.0 epoch
_J2000 = Time('J2000.0', scale='tt')


# One of the two nodes of intersection between the ICRF equator and Sun's equator in HCRS
_NODE = SkyCoord(_SOLAR_NORTH_POLE_HCRS.lon + 90*u.deg, 0*u.deg, frame='hcrs')


# The longitude in the de-tilted frame of the Sun's prime meridian.
# The IAU (Seidelmann et al. 2007 and later) defines the true longitude of the meridian (i.e.,
# without light travel time to Earth and aberration effects) as 84.176 degrees eastward at J2000.
_DLON_MERIDIAN = Longitude(_detilt_lon(_NODE) + constants.get('W_0'))


@add_common_docstring(**_variables_for_parse_time_docstring())
def L0(time='now',
        light_travel_time_correction=True,
        nearest_point=True,
        aberration_correction=False):
    """
    Return the L0 angle for the Sun at a specified time, which is the apparent Carrington longitude
    of the Sun-disk center as seen from Earth.

    Observer corrections can be disabled, and then this function will instead return the true
    Carrington longitude.

    Parameters
    ----------
    time : {parse_time_types}
        Time to use in a parse_time-compatible format
    light_travel_time_correction : `bool`
        If True, apply the correction for light travel time from Sun to Earth.  Defaults to True.
    nearest_point : `bool`
        If True, calculate the light travel time to the nearest point on the Sun's surface rather
        than the light travel time to the center of the Sun (i.e., a difference of the solar
        radius).  Defaults to True.
    aberration_correction : `bool`
        If True, apply the stellar-aberration correction due to Earth's motion.  Defaults to False.

    Returns
    -------
    `~astropy.coordinates.Longitude`
        The Carrington longitude

    Notes
    -----
    This longitude is calculated using current IAU values (Seidelmann et al. 2007 and later), which
    do not include the effects of light travel time and aberration due to Earth's motion (see that
    paper's Appendix).  This function then, by default, applies the light-travel-time correction
    for the nearest point on the Sun's surface, but does not apply the stellar-aberration correction
    due to Earth's motion.

    We do not apply the stellar-abberation correction by default because it should not be applied
    for purposes such as co-aligning images that are each referenced to Sun-disk center.  Stellar
    aberration does not shift the apparent positions of solar features relative to the Sun-disk
    center.

    The Astronomical Almanac applies the stellar-aberration correction in their printed published
    L0 values (see also Urban & Kaplan 2007).  Applying the stellar-aberration correction due to
    Earth's motion decreases the apparent Carrington longitude by ~20.5 arcseconds.

    References
    ----------
    * Seidelmann et al. (2007), "Report of the IAU/IAG Working Group on cartographic coordinates
      and rotational elements: 2006" `(link) <http://dx.doi.org/10.1007/s10569-007-9072-y>`__
    * Urban & Kaplan (2007), "Investigation of Change in the Computational Technique of the Sun's
      Physical Ephemeris in The Astronomical Almanac"
      `(link) <http://asa.hmnao.com/static/files/sun_rotation_change.pdf>`__
    """
    obstime = parse_time(time)
    earth = get_earth(obstime)

    # Calculate the de-tilt longitude of the Earth
    dlon_earth = _detilt_lon(earth)

    # Calculate the distance to the nearest point on the Sun's surface
    distance = earth.radius - constants.radius if nearest_point else earth.radius

    # Apply a correction for aberration due to Earth motion
    # This expression is an approximation to reduce computations (e.g., it does not account for the
    # inclination of the Sun's rotation axis relative to the ecliptic), but the estimated error is
    # <0.2 arcseconds
    if aberration_correction:
        dlon_earth -= 20.496*u.arcsec * 1*u.AU / earth.radius

    # Antedate the observation time to account for light travel time for the Sun-Earth distance
    antetime = (obstime - distance / speed_of_light) if light_travel_time_correction else obstime

    # Calculate the de-tilt longitude of the meridian due to the Sun's sidereal rotation
    dlon_meridian = Longitude(_DLON_MERIDIAN + (antetime - _J2000)
                              * constants.sidereal_rotation_rate)

    return Longitude(dlon_earth - dlon_meridian)


@add_common_docstring(**_variables_for_parse_time_docstring())
def P(time='now'):
    """
    Return the position (P) angle for the Sun at a specified time, which is the angle between
    geocentric north and solar north as seen from Earth, measured eastward from geocentric north.
    The range of P is +/-26.3 degrees.

    Parameters
    ----------
    time : {parse_time_types}
        Time to use in a parse_time-compatible format

    Returns
    -------
    out : `~astropy.coordinates.Angle`
        The position angle
    """
    obstime = parse_time(time)

    # Define the frame where its Z axis is aligned with geocentric north
    geocentric = ITRS(obstime=obstime)

    return _sun_north_angle_to_z(geocentric)


@add_common_docstring(**_variables_for_parse_time_docstring())
def earth_distance(time='now'):
    """
    Return the distance between the Sun and the Earth at a specified time.

    Parameters
    ----------
    time : {parse_time_types}
        Time to use in a parse_time-compatible format

    Returns
    -------
    out : `~astropy.coordinates.Distance`
        The Sun-Earth distance
    """
    obstime = parse_time(time)
    vector = get_body_barycentric('earth', obstime) - get_body_barycentric('sun', obstime)
    return Distance(vector.norm())


@add_common_docstring(**_variables_for_parse_time_docstring())
def orientation(location, time='now'):
    """
    Return the orientation angle for the Sun from a specified Earth location and time.  The
    orientation angle is the angle between local zenith and solar north, measured eastward from
    local zenith.

    Parameters
    ----------
    location : `~astropy.coordinates.EarthLocation`
        Observer location on Earth
    time : {parse_time_types}
        Time to use in a parse_time-compatible format

    Returns
    -------
    out : `~astropy.coordinates.Angle`
        The orientation of the Sun
    """
    obstime = parse_time(time)

    # Define the frame where its Z axis is aligned with local zenith
    local_frame = AltAz(obstime=obstime, location=location)

    return _sun_north_angle_to_z(local_frame)


def _sun_north_angle_to_z(frame):
    """
    Return the angle between solar north and the Z axis of the provided frame's coordinate system
    and observation time.
    """
    # Find the Sun center in HGS at the frame's observation time(s)
    sun_center_repr = SphericalRepresentation(0*u.deg, 0*u.deg, 0*u.km)
    # The representation is repeated for as many times as are in obstime prior to transformation
    sun_center = SkyCoord(sun_center_repr._apply('repeat', frame.obstime.size),
                          frame=HeliographicStonyhurst, obstime=frame.obstime)

    # Find the Sun north pole in HGS at the frame's observation time(s)
    sun_north_repr = SphericalRepresentation(0*u.deg, 90*u.deg, constants.radius)
    # The representation is repeated for as many times as are in obstime prior to transformation
    sun_north = SkyCoord(sun_north_repr._apply('repeat', frame.obstime.size),
                         frame=HeliographicStonyhurst, obstime=frame.obstime)

    # Find the Sun center and Sun north in the frame's coordinate system
    sky_normal = sun_center.transform_to(frame).data.to_cartesian()
    sun_north = sun_north.transform_to(frame).data.to_cartesian()

    # Use cross products to obtain the sky projections of the two vectors (rotated by 90 deg)
    sun_north_in_sky = sun_north.cross(sky_normal)
    z_in_sky = CartesianRepresentation(0, 0, 1).cross(sky_normal)

    # Normalize directional vectors
    sky_normal /= sky_normal.norm()
    sun_north_in_sky /= sun_north_in_sky.norm()
    z_in_sky /= z_in_sky.norm()

    # Calculate the signed angle between the two projected vectors
    cos_theta = sun_north_in_sky.dot(z_in_sky)
    sin_theta = sun_north_in_sky.cross(z_in_sky).dot(sky_normal)
    angle = np.arctan2(sin_theta, cos_theta).to('deg')

    # If there is only one time, this function's output should be scalar rather than array
    if angle.size == 1:
        angle = angle[0]

    return Angle(angle)
