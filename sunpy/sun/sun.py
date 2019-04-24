"""
This module provides Sun-related parameters.
"""
import numpy as np

import astropy.units as u
from astropy.coordinates import Angle, Latitude, Longitude

from sunpy.sun import constants
from sunpy.time import julian_centuries, parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring

__all__ = [
    "print_params", "apparent_declination",
    "apparent_rightascension", "apparent_obliquity_of_ecliptic", "true_declination",
    "true_rightascension", "true_obliquity_of_ecliptic", "apparent_latitude", "true_latitude",
    "apparent_longitude", "true_anomaly", "true_longitude",
    "equation_of_center", "geometric_mean_longitude", "carrington_rotation_number", "mean_anomaly",
    "longitude_sun_perigee", "mean_ecliptic_longitude", "eccentricity_sun_earth_orbit", "position",
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
    Return the angular size of the semi-diameter of the Sun as a function of
    time as viewed from Earth (in arcsec)

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
    celestial sphere using the equatorial coordinate system in arcsec.

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
def eccentricity_sun_earth_orbit(t='now'):
    """
    Returns the eccentricity of the Sun Earth Orbit.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a palongitude_Sun_perigee, true_latitude and apparent_latitude need to be fixedrse_time-compatible
        time string, number, or a datetime object.
    """
    T = julian_centuries(t)
    result = 0.016751040 - 0.00004180 * T - 0.0000001260 * T**2
    return result


@add_common_docstring(**_variables_for_parse_time_docstring())
def mean_ecliptic_longitude(t='now'):
    """
    Returns the mean ecliptic longitude.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    T = julian_centuries(t)
    result = 279.696680 + 36000.76892 * T + 0.0003025 * T**2
    result = result * u.deg
    return Longitude(result)


@add_common_docstring(**_variables_for_parse_time_docstring())
def longitude_sun_perigee(t='now'):
    """
    Returns the current solar perigee.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    # TODO: FIX THIS
    return 1


@add_common_docstring(**_variables_for_parse_time_docstring())
def mean_anomaly(t='now'):
    """
    Returns the mean anomaly (the angle through which the Sun has moved
    assuming a circular orbit) as a function of time.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    T = julian_centuries(t)
    result = 358.475830 + 35999.049750 * T - 0.0001500 * T**2 - 0.00000330 * T**3
    result = result * u.deg
    return Longitude(result)


@add_common_docstring(**_variables_for_parse_time_docstring())
def carrington_rotation_number(t='now'):
    """
    Return the Carrington Rotation number.

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
def geometric_mean_longitude(t='now'):
    """
    Returns the geometric mean longitude (in degrees).

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    T = julian_centuries(t)
    result = 279.696680 + 36000.76892 * T + 0.0003025 * T**2
    result = result * u.deg
    return Longitude(result)


@add_common_docstring(**_variables_for_parse_time_docstring())
def equation_of_center(t='now'):
    """
    Returns the Sun's equation of center (in degrees).

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    T = julian_centuries(t)
    mna = mean_anomaly(t)
    result = ((1.9194600 - 0.0047890 * T - 0.0000140 * T**2) * np.sin(mna) +
              (0.0200940 - 0.0001000 * T) * np.sin(2 * mna) + 0.0002930 * np.sin(3 * mna))
    result = result * u.deg
    return Angle(result)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_longitude(t='now'):
    """
    Returns the Sun's true geometric longitude (in degrees).

    Referred to the mean equinox of date.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    # TODO: Should the higher accuracy terms from which app_long is derived be added to true_long?
    result = equation_of_center(t) + geometric_mean_longitude(t)
    return Longitude(result)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_anomaly(t='now'):
    """
    Returns the Sun's true anomaly (in degrees).

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    result = mean_anomaly(t) + equation_of_center(t)
    return Longitude(result)


@add_common_docstring(**_variables_for_parse_time_docstring())
def apparent_longitude(t='now'):
    """
    Returns the apparent longitude of the Sun.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    T = julian_centuries(t)
    omega = (259.18 - 1934.142 * T) * u.deg
    true_long = true_longitude(t)
    result = true_long - (0.00569 - 0.00479 * np.sin(omega)) * u.deg
    return Longitude(result)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_latitude(t='now'):
    """
    Returns the true latitude.

    Never more than 1.2 arcsec from 0, set to 0 here.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    # TODO: FIX THIS
    return Latitude(0.0 * u.deg)


@add_common_docstring(**_variables_for_parse_time_docstring())
def apparent_latitude(t='now'):
    """
    Returns the true latitude.

    Set to 0 here.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    # TODO: FIX THIS
    return Latitude(0.0 * u.deg)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_obliquity_of_ecliptic(t='now'):
    """
    Returns the true obliquity of the ecliptic.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    T = julian_centuries(t)
    result = 23.452294 - 0.0130125 * T - 0.00000164 * T**2 + 0.000000503 * T**3
    return Angle(result, u.deg)


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_rightascension(t='now'):
    """
    Return the true right ascension.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    y = np.cos(true_obliquity_of_ecliptic(t)) * np.sin(true_longitude(t))
    x = np.cos(true_longitude(t))
    true_ra = np.arctan2(y, x)
    return Longitude(true_ra.to(u.hourangle))


@add_common_docstring(**_variables_for_parse_time_docstring())
def true_declination(t='now'):
    """
    Return the true declination.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    result = np.arcsin(np.sin(true_obliquity_of_ecliptic(t)) * np.sin(apparent_longitude(t)))
    return Latitude(result.to(u.deg))


@add_common_docstring(**_variables_for_parse_time_docstring())
def apparent_obliquity_of_ecliptic(t='now'):
    """
    Return the apparent obliquity of the ecliptic.

    Parameters
    ----------
    t : {parse_time_types}
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    """
    omega = apparent_longitude(t)
    result = true_obliquity_of_ecliptic(t) + (0.00256 * np.cos(omega)) * u.deg
    return result


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
    y = np.cos(apparent_obliquity_of_ecliptic(t)) * np.sin(apparent_longitude(t))
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
    ob = apparent_obliquity_of_ecliptic(t)
    app_long = apparent_longitude(t)
    result = np.arcsin(np.sin(ob)) * np.sin(app_long)
    return Latitude(result.to(u.deg))


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

    print('Solar Ephemeris for {}\n'.format(parse_time(t).ctime()))
    print('Distance = {}'.format(get_sunearth_distance(t)))
    print('Semidiameter = {}'.format(solar_semidiameter_angular_size(t)))
    print('True (long, lat) = ({}, {})'.format(true_longitude(t), true_latitude(t)))
    print('Apparent (long, lat) = ({}, {})'.format(apparent_longitude(t), apparent_latitude(t)))
    print('True (RA, Dec) = ({}, {})'.format(true_rightascension(t), true_declination(t)))
    print('Apparent (RA, Dec) = ({}, {})'.format(apparent_rightascension(t),
                                                 apparent_declination(t)))
    print('Heliographic long. and lat of disk center = ({}, {})'.format(get_sun_L0(t),
                                                                        get_sun_B0(t)))
    print('Position angle of north pole in = {}'.format(get_sun_P(t)))
    print('Carrington Rotation Number = {}'.format(carrington_rotation_number(t)))
