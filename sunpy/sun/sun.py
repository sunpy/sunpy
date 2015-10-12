"""Provides Sun-related parameters

The following code is heavily based on IDL function get_sun.pro which itself
is based on algorithms presented in the book Astronomical Formulae for
Calculators, by Jean Meeus.
Every function returning a quantity is of type astropy.units.Quantity

A correct answer set to compare to

Solar Ephemeris for  1-JAN-01  00:00:00

Distance (AU) = 0.98330468
Semidiameter (arc sec) = 975.92336
True (long, lat) in degrees = (280.64366, 0.00000)
Apparent (long, lat) in degrees = (280.63336, 0.00000)
True (RA, Dec) in hrs, deg = (18.771741, -23.012449)
Apparent (RA, Dec) in hrs, deg = (18.770994, -23.012593)
Heliographic long. and lat. of disk center in deg = (217.31269, -3.0416292)
Position angle of north pole in deg = 2.0102649
Carrington Rotation Number = 1971.4091        check!

"""
from __future__ import absolute_import, division, print_function

import numpy as np

import astropy.units as u
from astropy.coordinates import Angle, Longitude, Latitude

from sunpy.time import parse_time, julian_day, julian_centuries
from sunpy.sun import constants

__all__ = ["print_params"
           ,"heliographic_solar_center"
           ,"solar_north"
           ,"apparent_declination"
           ,"apparent_rightascension"
           ,"apparent_obliquity_of_ecliptic"
           ,"true_declination"
           ,"true_rightascension"
           ,"true_obliquity_of_ecliptic"
           ,"apparent_latitude"
           ,"true_latitude"
           ,"apparent_longitude"
           ,"sunearth_distance"
           ,"true_anomaly"
           ,"true_longitude"
           ,"equation_of_center"
           ,"geometric_mean_longitude"
           ,"carrington_rotation_number"
           ,"mean_anomaly"
           ,"longitude_Sun_perigee"
           ,"mean_ecliptic_longitude"
           ,"eccentricity_SunEarth_orbit"
           ,"position"
           ,"solar_semidiameter_angular_size"
           ,"solar_cycle_number"]

__authors__ = ["Steven Christe"]
__email__ = "steven.d.christe@nasa.gov"

def solar_cycle_number(t='now'):
    """Return the solar cycle number."""
    time = parse_time(t)
    result = (time.year + 8) % 28 + 1
    return result

def solar_semidiameter_angular_size(t='now'):
    r"""
    Return the angular size of the semi-diameter of the Sun as
    a function of time as viewed from Earth (in arcsec)

    .. math::

        Radius_{\odot}[rad]=\tan^{-1}\left(\frac{<Radius_{\odot}[m]>}{D_{\odot \oplus}(t)[m]}\right)

    since :math:`tan(x) \approx x` when :math:`x << 1`

    .. math::

        Radius_{\odot}[rad]=\frac{<Radius_{\odot}[m]>}{D_{\odot \oplus}(t)[m]}

    """
    solar_semidiameter_rad = (constants.radius.to(u.AU)) / sunearth_distance(t)
    return Angle(solar_semidiameter_rad.to(u.arcsec, equivalencies = u.dimensionless_angles()))

def position(t='now'):
    """Returns the position of the Sun (right ascension and declination)
    on the celestial sphere using the equatorial coordinate system in arcsec.
    """
    ra = true_rightascension(t)
    dec = true_declination(t)
    result = [ra,dec]
    return result

def eccentricity_SunEarth_orbit(t='now'):
    """Returns the eccentricity of the Sun Earth Orbit."""
    T = julian_centuries(t)
    result = 0.016751040 - 0.00004180 * T - 0.0000001260 * T ** 2
    return result

def mean_ecliptic_longitude(t='now'):
    """Returns the mean ecliptic longitude."""
    T = julian_centuries(t)
    result = 279.696680 + 36000.76892 * T + 0.0003025 * T ** 2
    result = result * u.deg
    return Longitude(result)

def longitude_Sun_perigee(t='now'): # pylint: disable=W0613
    # T = julian_centuries(t)
    return 1

def mean_anomaly(t='now'):
    """Returns the mean anomaly (the angle through which the Sun has moved
    assuming a circular orbit) as a function of time."""
    T = julian_centuries(t)
    result = 358.475830 + 35999.049750 * T - 0.0001500 * T ** 2 - 0.00000330 * T ** 3
    result = result * u.deg
    return Longitude(result)

def carrington_rotation_number(t='now'):
    """Return the Carrington Rotation number"""
    jd = julian_day(t)
    result = (1. / 27.2753) * (jd - 2398167.0) + 1.0
    return result

def geometric_mean_longitude(t='now'):
    """Returns the geometric mean longitude (in degrees)"""
    T = julian_centuries(t)
    result = 279.696680 + 36000.76892 * T + 0.0003025 * T ** 2
    result = result * u.deg
    return Longitude(result)

def equation_of_center(t='now'):
    """Returns the Sun's equation of center (in degrees)"""
    T = julian_centuries(t)
    mna = mean_anomaly(t)
    result = ((1.9194600 - 0.0047890 * T - 0.0000140 * T ** 2) * np.sin(mna)
    + (0.0200940 - 0.0001000 * T) *
    np.sin(2 * mna) + 0.0002930 * np.sin(3 * mna))
    result = result * u.deg
    return Angle(result)

def true_longitude(t='now'):
    """Returns the Sun's true geometric longitude (in degrees)
    (Referred to the mean equinox of date.  Question: Should the higher
    accuracy terms from which app_long is derived be added to true_long?)"""
    result = equation_of_center(t) + geometric_mean_longitude(t)
    return Longitude(result)

def true_anomaly(t='now'):
    """Returns the Sun's true anomaly (in degrees)."""
    result = (mean_anomaly(t) + equation_of_center(t)) % (360.0 * u.deg)
    return result

def sunearth_distance(t='now'):
    """Returns the Sun Earth distance (AU). There are a set of higher
    accuracy terms not included here."""
    ta = true_anomaly(t)
    e = eccentricity_SunEarth_orbit(t)
    result = 1.00000020 * (1.0 - e ** 2) / (1.0 + e * np.cos(ta))
    return result * u.AU

def apparent_longitude(t='now'):
    """Returns the apparent longitude of the Sun."""
    T = julian_centuries(t)
    omega = (259.18 - 1934.142 * T) * u.deg
    true_long = true_longitude(t)
    result = true_long - (0.00569 - 0.00479 * np.sin(omega)) * u.deg
    return Longitude(result)

def true_latitude(t='now'): # pylint: disable=W0613
    """Returns the true latitude. Never more than 1.2 arcsec from 0,
    set to 0 here."""
    return 0.0

def apparent_latitude(t='now'): # pylint: disable=W0613
    """Returns the true latitude. Set to 0 here."""
    return 0

def true_obliquity_of_ecliptic(t='now'):
    """Returns the true obliquity of the ecliptic."""
    T = julian_centuries(t)
    result = 23.452294 - 0.0130125 * T - 0.00000164 * T ** 2 + 0.000000503 * T ** 3
    return Angle(result, u.deg)

def true_rightascension(t='now'):
    """Return the true right ascension."""
    true_long = true_longitude(t)
    ob = true_obliquity_of_ecliptic(t)
    result = np.cos(ob) * np.sin(true_long)
    result = result * u.deg
    return Longitude(result)

def true_declination(t='now'):
    """Return the true declination."""
    result = np.cos(true_longitude(t))
    result = result * u.deg
    return Latitude(result)

def apparent_obliquity_of_ecliptic(t='now'):
    """Return the apparent obliquity of the ecliptic."""
    omega = apparent_longitude(t)
    result = true_obliquity_of_ecliptic(t) + (0.00256 * np.cos(omega)) * u.deg
    return result

def apparent_rightascension(t='now'):
    """Returns the apparent right ascension of the Sun."""
    y = np.cos(apparent_obliquity_of_ecliptic(t)) * np.sin(apparent_longitude(t))
    x = np.cos(apparent_longitude(t))
    app_ra = np.arctan2(y, x)
    return Longitude(app_ra.to(u.hourangle))

def apparent_declination(t='now'):
    """Returns the apparent declination of the Sun."""
    ob = apparent_obliquity_of_ecliptic(t)
    app_long = apparent_longitude(t)
    result = np.degrees(np.arcsin(np.sin(ob)) * np.sin(app_long))
    return Latitude(result)

def solar_north(t='now'):
    """Returns the position of the Solar north pole in degrees."""
    T = julian_centuries(t)
    ob1 = true_obliquity_of_ecliptic(t)
    # in degrees
    i = 7.25 * u.deg
    k = (74.3646 + 1.395833 * T) * u.deg
    lamda = true_longitude(t) - (0.00569 * u.deg)
    omega = apparent_longitude(t)
    lamda2 = lamda - (0.00479 * np.sin(omega)) * u.deg
    diff = lamda - k
    x = np.arctan(-np.cos((lamda2) * np.tan(ob1)))
    y = np.arctan(-np.cos(diff) * np.tan(i))
    result = x + y
    return Angle(result.to(u.deg))

def heliographic_solar_center(t='now'):
    """Returns the position of the solar center in heliographic coordinates."""
    jd = julian_day(t)
    T = julian_centuries(t)
    # Heliographic coordinates in degrees
    theta = ((jd - 2398220)*360/25.38) * u.deg
    i = 7.25 * u.deg
    k = (74.3646 + 1.395833 * T) * u.deg
    lamda = true_longitude(t) - 0.00569 * u.deg
    diff = lamda - k
    # Latitude at center of disk (deg):
    he_lat = np.degrees(np.arcsin(np.sin(diff)*np.sin(i)))
    # Longitude at center of disk (deg):
    y = -np.sin(diff)*np.cos(i)
    x = -np.cos(diff)
    rpol = (np.arctan2(y, x))
    he_lon = rpol - theta
    return [Longitude(he_lon), Latitude(he_lat)]

def print_params(t='now'):
    """Print out a summary of Solar ephemeris"""
    time = parse_time(t)
    print('Solar Ephemeris for ' + time.ctime())
    print('')
    print('Distance (AU) = ' + str(sunearth_distance(t)))
    print('Semidiameter (arc sec) = ' + str(solar_semidiameter_angular_size(t)))
    print('True (long,lat) in degrees = (' + str(true_longitude(t)) + ','
                                                 + str(true_latitude(t)) + ')')
    print('Apparent (long, lat) in degrees = (' + str(apparent_longitude(t)) + ','
                                                 + str(apparent_latitude(t)) + ')')
    print('True (RA, Dec) = (' + str(true_rightascension(t)) + ','
          + str(true_declination(t)))
    print('Apparent (RA, Dec) = (' + str(apparent_rightascension(t)) + ','
          + str(apparent_declination(t)))
    print('Heliographic long. and lat of disk center in deg = (' + str(heliographic_solar_center(t)) + ')')
    print('Position angle of north pole in deg = ' + str(solar_north(t)))
    print('Carrington Rotation Number = ' + str(carrington_rotation_number(t)))
