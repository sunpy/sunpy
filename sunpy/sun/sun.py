"""
This module provides Sun-related parameters.
"""
import numpy as np

import astropy.units as u
from astropy.coordinates import Angle, Latitude, Longitude, SkyCoord
# Versions of Astropy that do not have GeocentricMeanEcliptic have the same frame
# with the incorrect name GeocentricTrueEcliptic
try:
    from astropy.coordinates import GeocentricMeanEcliptic
except ImportError:
    from astropy.coordinates import GeocentricTrueEcliptic as GeocentricMeanEcliptic
from astropy import _erfa as erfa
from astropy.coordinates.builtin_frames.utils import get_jd12

from sunpy.sun import constants
from sunpy.time import julian_centuries, parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring

__all__ = [
    "print_params", "apparent_declination",
    "apparent_rightascension", "true_obliquity_of_ecliptic", "true_declination",
    "true_rightascension", "mean_obliquity_of_ecliptic", "apparent_latitude", "true_latitude",
    "apparent_longitude", "true_longitude", "carrington_rotation_number", "position",
    "solar_semidiameter_angular_size", "solar_cycle_number"
]


@add_common_docstring(**_variables_for_parse_time_docstring())
def solar_cycle_number(t='now'):
    """
    Return the solar cycle number.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    time = parse_time(t)
    result = (int(time.strftime('%Y')) + 8) % 28 + 1
    return result


@add_common_docstring(**_variables_for_parse_time_docstring())
def solar_semidiameter_angular_size(t='now'):
    """
    Return the angular size of the semi-diameter of the Sun as viewed from Earth.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    # Import here to avoid a circular import
    from sunpy.coordinates import get_sunearth_distance
    solar_semidiameter_rad = (constants.radius.to(u.AU)) / get_sunearth_distance(t)
    return Angle(solar_semidiameter_rad.to(u.arcsec, equivalencies=u.dimensionless_angles()))


@add_common_docstring(**_variables_for_parse_time_docstring())
def position(t='now'):
    """
    Returns the position of the Sun (right ascension and declination) on the
    celestial sphere using the equatorial coordinate system.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    ra = true_rightascension(t)
    dec = true_declination(t)
    return ra, dec


@add_common_docstring(**_variables_for_parse_time_docstring())
def carrington_rotation_number(t='now'):
    """
    Return the Carrington Rotation number.

    .. warning::
        The accuracy of this function is under investigation.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    jd = parse_time(t).jd
    result = (1. / 27.2753) * (jd - 2398167.0) + 1.0
    return result


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_longitude(t='now'):
    """
    Returns the Sun's true geometric longitude, referred to the mean equinox of date.  No
    corrections for nutation or aberration are included.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    time = parse_time(t)
    sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='hcrs', obstime=time)
    coord = sun.transform_to(GeocentricMeanEcliptic(equinox=time))
    # Astropy's GeocentricMeanEcliptic includes aberration from Earth motion, so remove it
    lon = coord.lon + Angle('20.496s') * 1*u.AU / coord.distance
    return Longitude(lon)


@add_common_docstring(**_variables_for_parse_time_docstring())
def apparent_longitude(t='now'):
    """
    Returns the Sun's apparent longitude, referred to the true equinox of date.  Corrections for
    nutation and aberration (for Earth motion) are included.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    time = parse_time(t)
    sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='hcrs', obstime=time)
    coord = sun.transform_to(GeocentricMeanEcliptic(equinox=time))

    # Astropy's GeocentricMeanEcliptic already includes aberration, so only add nutation
    jd1, jd2 = get_jd12(time, 'tt')
    nut_lon, _ = erfa.nut06a(jd1, jd2)*u.radian
    lon = coord.lon + nut_lon

    return Longitude(lon)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_latitude(t='now'):
    """
    Returns the Sun's true geometric latitude, referred to the mean equinox of date.  No
    corrections for nutation or aberration are included.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    time = parse_time(t)
    sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='hcrs', obstime=time)
    coord = sun.transform_to(GeocentricMeanEcliptic(equinox=time))

    # Astropy's GeocentricMeanEcliptic includes aberration from Earth motion, but it's ignorable
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
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    time = parse_time(t)
    sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame='hcrs', obstime=time)
    coord = sun.transform_to(GeocentricMeanEcliptic(equinox=time))

    # Astropy's GeocentricMeanEcliptic does not include nutation, but it's ignorable
    lat = coord.lat

    return Latitude(lat)


@add_common_docstring(**_variables_for_parse_time_docstring())
def mean_obliquity_of_ecliptic(t='now'):
    """
    Returns the mean obliquity of the ecliptic, using the IAU 2006 definition.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    time = parse_time(t)
    jd1, jd2 = get_jd12(time, 'tt')
    obl = erfa.obl06(jd1, jd2)*u.radian
    return Angle(obl, u.arcsec)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_rightascension(t='now'):
    """
    Return the true right ascension of the Sun.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    y = np.cos(mean_obliquity_of_ecliptic(t)) * np.sin(true_longitude(t))
    x = np.cos(true_longitude(t))
    true_ra = np.arctan2(y, x)
    return Longitude(true_ra, u.hourangle)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_declination(t='now'):
    """
    Return the true declination of the Sun.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    result = np.arcsin(np.sin(mean_obliquity_of_ecliptic(t)) * np.sin(true_longitude(t)))
    return Latitude(result, u.deg)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_obliquity_of_ecliptic(t='now'):
    """
    Returns the true obliquity of the ecliptic, using the IAU 2006 definition.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    time = parse_time(t)
    jd1, jd2 = get_jd12(time, 'tt')
    obl = erfa.obl06(jd1, jd2)*u.radian
    _, nut_obl = erfa.nut06a(jd1, jd2)*u.radian
    obl += nut_obl
    return Angle(obl, u.arcsec)


@add_common_docstring(**_variables_for_parse_time_docstring())
def apparent_rightascension(t='now'):
    """
    Returns the apparent right ascension of the Sun.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    y = np.cos(true_obliquity_of_ecliptic(t)) * np.sin(apparent_longitude(t))
    x = np.cos(apparent_longitude(t))
    app_ra = np.arctan2(y, x)
    return Longitude(app_ra.to(u.hourangle))


@add_common_docstring(**_variables_for_parse_time_docstring())
def apparent_declination(t='now'):
    """
    Returns the apparent declination of the Sun.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    result = np.arcsin(np.sin(true_obliquity_of_ecliptic(t)) * np.sin(apparent_longitude(t)))
    return Latitude(result, u.deg)


@add_common_docstring(**_variables_for_parse_time_docstring())
def print_params(t='now'):
    """
    Print out a summary of Solar ephemeris.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    # import here to avoid circular import
    from sunpy.coordinates.ephemeris import (get_sun_L0, get_sun_B0,
                                             get_sun_P, get_sunearth_distance)

    print('Solar Ephemeris for {} UTC\n'.format(parse_time(t).utc))
    print('Distance = {}'.format(get_sunearth_distance(t)))
    print('Semidiameter = {}'.format(solar_semidiameter_angular_size(t)))
    print('True (long, lat) = ({}, {})'.format(true_longitude(t).to_string(),
                                               true_latitude(t).to_string()))
    print('Apparent (long, lat) = ({}, {})'.format(apparent_longitude(t).to_string(),
                                                   apparent_latitude(t).to_string()))
    print('True (RA, Dec) = ({}, {})'.format(true_rightascension(t).to_string(),
                                             true_declination(t).to_string()))
    print('Apparent (RA, Dec) = ({}, {})'.format(apparent_rightascension(t).to_string(),
                                                 apparent_declination(t).to_string()))
    print('Heliographic long. and lat of disk center = ({}, {})'.format(get_sun_L0(t).to_string(),
                                                                        get_sun_B0(t).to_string()))
    print('Position angle of north pole in = {}'.format(get_sun_P(t)))
    print('Carrington Rotation Number = {}'.format(carrington_rotation_number(t)))
