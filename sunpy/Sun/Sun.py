# -*- coding: utf-8 -*-
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
#
# <License info will go here...>
"""Provides Sun-related parameters

    Following code is heavily based on IDL function get_sun.pro
    which itself is based on algorithms presented in the book
    Practical Astronomy with your Calculator (Third Edition)
    by Peter Duffet-Smith

    Notes:

    Keith (2011/03/29)
     If possible, it might be helpful to encode some of the magic
     numbers as constant, e.g. CONSTANT_NAME = 12345.6     
"""

import datetime
import math
import cmath
import numpy as np
#TODO: why do i need to use util.util?
import sunpy.util.util as util

def radius(t=None):
    return angular_size(t)

def angular_size(t=None):
    """Return the angular size of the Sun as a function of 
    time as viewed from Earth (in arcsec)"""
    result = 959.63 / sunearth_distance(t)
    return result

def position(t=None):
	"""Returns the position of the Sun (right ascension and declination)
	on the celestial sphere using the equatorial coordinate system in arcsec.
	"""
	return 1
    
def eccentricity_SunEarth_orbit(t=None):
	"""Returns the eccentricity of the Sun Earth Orbit."""
	T = util.julian_centuries(t)
	result = 0.016751040 - 0.00004180 * T - 0.0000001260 * T ** 2
	return result

def mean_ecliptic_longitude(t=None):
    """Returns the mean ecliptic longitude."""
    T = util.julian_centuries(t)
    result = 279.696680 + 36000.76892 * T + 0.0003025 * T ** 2
    result = result % 360.0
    return result

def longitude_Sun_perigee(t=None):
    """Insert text here"""
    T = util.julian_centuries(t)
    return 1
	
def mean_anomaly(t=None):
    """Returns the mean anomaly (the angle through which the Sun has moved
    assuming a circular orbit) as a function of time."""
    T = util.julian_centuries(t)
    result = 358.475830 + 35999.049750 * T - 0.0001500 * T ** 2 - 0.00000330 * T ** 3
    result = result % 360.0
    return result

def carrington_rotation_number(t=None):
    """Return the Carrington Rotation number"""
    jd = julian_day(t)
    result = (1. / 27.2753) * (jd - 2398167.0) + 1.0
    return result

def geometric_mean_longitude(t=None):
    """Returns the geometric mean longitude (in degrees)"""   
    T = util.julian_centuries(t)
    result = 279.696680 + 36000.76892 * T + 0.0003025 * T ** 2
    result = result % 360.0
    return result
  
def equation_of_center(t=None):
    """Returns the Sun's equation of center (in degrees)"""
    T = util.julian_centuries(t)
    mna = mean_anomaly(t) 
    result = ((1.9194600 - 0.0047890 * T - 0.0000140 * T
    ** 2) * np.sin(np.radians(mna) + (0.0200940 - 0.0001000 * T) *
    np.sin(np.radians(2 * mna)) + 0.0002930 * np.sin(np.radians(3 * mna))))
    return result

def geometric_longitude(t=None): 
    """Returns the Sun's true geometric longitude (in degrees) 
    (Refered to the mean equinox of date.  Question: Should the higher
    accuracy terms from which app_long is derived be added to true_long?)"""
    result = (equation_of_center(util.anytitomaly(t))) % 360.0
    return result

def true_anomaly(t=None):
    """Returns the Sun's true anomaly (in degress)."""
    result = (mean_anomaly(t) + equation_of_center(t)) % 360.0
    return result

def sunearth_distance(t=None):
    """Returns the Sun Earth distance. There are a set of higher accuracy terms not included here."""  
    ta = true_anomaly(t)
    e = eccentricity_SunEarth_orbit(t)
    result = 1.00000020 * (1.0 - e ** 2) / (1.0 + e * math.cos(np.radians(ta)))
    return result

def apparent_longitude(t=None):
    """Returns the apparent longitude of the Sun."""
    T = util.julian_centuries(t)
    omega = 259.18 - 1934.142*T
    true_long = geometric_longitude(t)        
    result = true_long - 0.00569 - 0.00479*math.sin(np.radian(omega))
    return result

def geometric_latitude(t=None):
    return 0

def apparent_latitude(t=None):
    return 0

def true_obliquity_of_ecliptic(t=None):
    T = util.julian_centuries(t)
    result = 23.452294 - 0.0130125 * T - 0.00000164 * T ** 2 + 0.000000503 * t ** 3
    return result

def true_rightascenscion(t=None):
    true_long = geometric_longitude(t)
    ob = true_obliquity_of_ecliptic(t)
    result = math.cos(np.radian(ob))*math.sin(np.radian(true_long))
    return result

def true_declination(t=None):
    result = math.cos(np.radian(geometric_longitude(t)))
    return result

def apparent_obliquity_of_ecliptic(t=None):
    result = true_obliquity_of_ecliptic(t) + 0.00256 * math.cos(np.radian(omega))
    return result

def apparent_rightascenscion(t=None):
    """Returns the apparent right ascenscion of the Sun."""
    return result

def apparent_declination(t=None):
    """Returns the apparent declination of the Sun."""
    ob = apparent_obliquity_of_ecliptic(t)
    app_long = apparent_longitude(t)
    result = np.degrees(math.asin(math.sin(np.radian(ob)))*math.sin(np.radian(app_long)))
    return result

def solar_north(t=None):
    """Returns the position of the Solar north pole in degrees."""
    jd = julian_day(t)
    T = util.julian_centuries(t)

    theta = (jd - 2398220)*360/25.38
    # in degrees
    i = 7.25
    k = 74.3646 + 1.395833 * T
    lamda = true_longitude(t) - 0.00569
    lamda2 = lamda - 0.00479 * math.sin(np.radian(omega))
    diff = np.radian(lamda - k)
    x = np.degrees(math.atan(-math.cos(np.radian(lamda2)*math.tan(np.radian(ob1)))))
    y = np.degrees(math.atan(-math.cos(diff)*math.tan(np.radian(i))))
    result = x + y
    return result

def heliographic_solar_center(t=None):
    """Returns the position of the solar center in heliographic coordinates."""
    jd = julian_day(t)
    T = util.julian_centuries(t)
    # Heliographic coordinates in degrees
    theta = (jd - 2398220)*360/25.38
    i = 7.25
    k = 74.3646 + 1.395833 * T
    lamda = true_longitude(t) - 0.00569
    lamda2 = lamda - 0.00479 * math.sin(np.radian(omega))
    diff = np.radian(lamda - k)
    # Latitude at center of disk (deg):    
    he_lat = np.degrees(math.asin(math.sin(diff)*sin(np.radian(i))))
    # Longitude at center of disk (deg):
    y = -math.sin(diff)*cos(np.radian(i))
    x = -math.cos(diff)
    rpol = cmath.polar(complex(x,y))
    he_lon = np.degrees(rpol[1]) - theta
    he_lon = he_lon % 360
    if he_lon < 0:
        he_lon = he_lon + 360.0

    return [he_lon, he_lat]

def test(t=None):
    print(angular_size(t))
    
