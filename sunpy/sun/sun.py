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
from __future__ import absolute_import

import numpy as np

import astropy.units as u
import astropy.time
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
    """
    Returns approximately the solar cycle number from 1755 assuming an average
    of 10.66 years per cycle. The average has been obtained having into account
    Usoskin et al. (2009) results about the total number of cycles between
    1699 till 2008.

    Usoskin et al. (2009) -  doi:10.1088/0004-637X/700/2/L154
    """
    time = parse_time(t)
    result = int((time.year - 1755) * u.yr / constants.constant('sunspot cycle'))
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
    """Returns the eccentricity of the Sun Earth Orbit.
    Eq. 25.4
    """
    #TODO: Is this low accuracy?
    T = julian_centuries(t)
    result = 0.016708634 - 0.000042037 * T - 0.0000001267 * T ** 2
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
    assuming a circular orbit) as a function of time referred to the mean
    equinox of the date.
    The mean anomaly of the Sun is the same as the mean anomaly of the Earth
    (Meeus, 2005)
    """
    T = julian_centuries(t)
    mean_longitude_earth = np.array([100.466457, 36000.7698278, 3.0322e-4,
                                     +2.0e-8]) * u.deg
    lon_perihelium_earth = np.array([102.937348,     1.7195366, 4.5688e-4,
                                     -1.8e-8]) * u.deg
    polinomial_order = np.array([[1, T, T**2, T**3]])
    mean_anomaly = (mean_longitude_earth - lon_perihelium_earth) * polinomial_order
    return Longitude(mean_anomaly.sum())

def carrington_rotation_number(t='now'):
    """Return the Carrington Rotation number"""
    jd = julian_day(t)
    result = (1. / 27.2753) * (jd - 2398167.0) + 1.0
    return result

def geometric_mean_longitude(t='now'):
    """Returns the geometric mean longitude (in degrees)"""   
    T = julian_centuries(t)
    #Eq 25.2
    gm_long = np.array([280.46646, 36000.76983, 0.0003032]) * u.deg
    gm_long *= T ** np.arange(3)
    #TODO: High accuracy?
    return Longitude(gm_long)
  
def equation_of_center(t='now'):
    """Returns the Sun's equation of center (in degrees)
    p.164 """
    T = julian_centuries(t)
    mna = mean_anomaly(t) 
    factors = np.array([1.914602 - 0.004817 * T - 0.0000140 * T ** 2,
                        0.019993 - 0.000101 * T,
                        0.000289]) * u.deg
    result = (factors * np.sin(np.arange(1,4) * mna)).sum()
    result = result * u.deg
    #TODO: is it highaccuracy?
    return Angle(result)

def true_longitude(t='now'):
    """Returns the Sun's true geometric longitude (in degrees) 
    (Refered to the mean equinox of date.  Question: Should the higher
    accuracy terms from which app_long is derived be added to true_long?)"""
    #p.164
    result = geometric_mean_longitude(t) + equation_of_center(t)
    # Should we use p166 (VSOP87 or FK5)
    return Longitude(result) #TODO: print gives it in days, minutes, seconds

def true_anomaly(t='now'):
    """Returns the Sun's true anomaly (in degress)."""
    #p.164
    result = mean_anomaly(t) + equation_of_center(t)
    return Angle(result)

def sunearth_distance(t='now'):
    """Returns the Sun Earth distance (AU).
    #TODO: There are a set of higher 
    accuracy terms not included here.
    Eq 25.5 """  
    ta = true_anomaly(t)
    e = eccentricity_SunEarth_orbit(t)
    result = 1.000001018 * (1.0 - e ** 2) / (1.0 + e * np.cos(ta))
    return result * u.AU

def apparent_longitude(t='now', high_precission=False):
    """Returns the apparent longitude of the Sun.
    This includes the effects of nutation and aberration 
    to the true ('geometric') longituded of the Sun."""
    # p. 167
    # * Nutation: Add \sun longitude the nutation in long \delta\psi (Ch22)
    # ToDo: find!
    # * Aberration: Apply to \sun longitude the correction:
    # - 20.4898 * u.arcsecs / R
    # The numerator = constant of aberration (k=20.49552 *u.arcsecs)
    #                 * a(1-e**2) = eq. 25.10
    # It varies very slowly
    # Note, Above: it does not assume perturbations in the orbit (mainly by
    #       the moon => result's error <= 0.01 arcseconds
    # --High accuracy--
    # * Aberration: 
    # - 0.005775518 R \delta\lambda (25.11)
    #   the numerical const: light-time for unit distance, in days (=8.3 min).
    #   \deltlambda (J2000.0) = sun_geocentric_daily_variation
    # most important periodic terms have been retained => err < 0.1 arcsecs
    #  If this result is used for 25.11 => err < 0.001 arcsecs
    
    T = julian_centuries(t)
    true_long = true_longitude(t)
    if high_precission:
        nutation = 0 #ToDo CH22
        R = sunearth_distance(t)
        aberration = -0.005775518 * R * sun_geocentric_daily_variation(t)
        result = true_long + nutation + aberration
    else:
        omega = (259.18 - 1934.142 * T) * u.deg
        result = true_long - (0.00569 - 0.00479 * np.sin(omega)) * u.deg
    return Longitude(result)

def sun_geocentric_daily_variation(t='now', mean_equinox=False):
    """
    Daily variation, in arcseconds, of the geocentric longitude
    of the Sun in a fixed reference frame.
    If needed with respect to the mean equinox instead of to a 
    fixed reference frame, set argument to True.
    """
    # tau is measured from J2000.0 (JDE 2451 545.0) in Julian millennia
    tau = julian_centuries(t)/10.

    #  If \deltaLambda needed with respect to the mean equinox of the
    # date instead of to a fixed reference frame the constant term
    # should be replaced: 3548.193 -> 3548.330
    # p.167
    if mean_equinox:
        deltaLambda = 3548.330
    else:
        deltaLambda = 3548.193

    # Periodic terms due to Earth's orbit eccentricity
    earth_pt = np.array([359993.7286, 719987.4571, 1079981.1857])
    #                due to the action of the Moon
    moon_pt = np.array([4452671.1152])
    #                due to Venus
    venus_pt = np.array([450368.8564, 225184.4282, 315559.5560, 675553.2846])
    #                due to Jupiter
    jupit_pt = np.array([329644.6718, 659289.3436, 299295.6151])
    #                due to Mars
    mars_pt = np.array([337181.4711])
    periodic_terms = np.array([ earth_pt[0],  earth_pt[1], moon_pt[0],
                                venus_pt[0],  jupit_pt[0], jupit_pt[1],
                                9224659.7915, earth_pt[2], venus_pt[1],
                                4092677.3866, mars_pt[0],  jupit_pt[2],
                                venus_pt[2],  venus_pt[3], earth_pt[0],
                                earth_pt[1],  earth_pt[2], earth_pt[0],
                                earth_pt[1],  moon_pt[0],  earth_pt[0]])
    non_periodic_terms = np.array([  87.5287,  85.0561,  27.8502,
                                     73.1375, 337.2264, 222.5400,
                                    162.8136,  82.5823, 171.5189,
                                     30.3214, 119.8105, 247.5418,
                                    325.1526, 155.1241, 333.4515,
                                    330.9814, 328.5170, 241.4518,
                                    205.0482, 297.8610, 154.7066])
    factor = np.array([ 118.568, 2.476, 1.376,
                          0.119, 0.114, 0.086,
                          0.078, 0.054, 0.052,
                          0.034, 0.033, 0.023,
                          0.023, 0.021,
                          7.311, 0.305, 0.010, # * tau
                          0.309, 0.021, 0.004, # * tau ** 2
                          0.010])              # * tau ** 3
    factor[14:17] *= tau
    factor[17:20] *= tau ** 2
    factor[-1] *= tau ** 3

    arguments = factor * np.sin((non_periodic_terms + periodic_terms * tau) *
                                u.deg)
    deltaLambda += arguments.sum()
    return deltaLambda * u.arcsec 

def true_latitude(t='now'): # pylint: disable=W0613
    '''Returns the true latitude. Never more than 1.2 arcsec from 0,
    set to 0 here.'''
    return 0.0

def apparent_latitude(t='now'): # pylint: disable=W0613
    return 0

def true_obliquity_of_ecliptic(t='now'):
    T = julian_centuries(t)
    result = 23.452294 - 0.0130125 * T - 0.00000164 * T ** 2 + 0.000000503 * T ** 3
    return Angle(result, u.deg)

def true_rightascension(t='now'):
    true_long = true_longitude(t)
    ob = true_obliquity_of_ecliptic(t)
    result = np.cos(ob) * np.sin(true_long)
    result = result * u.deg
    return Longitude(result)

def true_declination(t='now'):
    result = np.cos(true_longitude(t))
    result = result * u.deg
    return Latitude(result)

def apparent_obliquity_of_ecliptic(t='now'):
    omega = apparent_longitude(t)
    result = true_obliquity_of_ecliptic(t) + (0.00256 * np.cos(omega)) * u.deg
    return result

def apparent_rightascension(t='now'):
    """Returns the apparent right ascenscion of the Sun."""
    # P.167: apparent long (\lambda) and lat (\beta) can be transformed
    # to app_ra and app_dec (13.3) and (13.4) where \eps = true obliquity
    # affected by nutation in obliquitiy (\delta \eps)
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

def longitude_ascending_node(t='now'):
    """Returns the longitude of the ascending node of the solar equator on
    the ecliptic.

    Astronomical Algorithms 2nd Ed. - Jean Meeus 2005  - ISBN: 0-943396-61-1
    """
    ut = astropy.time.Time(parse_time(t), scale = 'utc')
    jd_ephem = ut.tt.jd
    k = 73.6667 * u.deg + 1.3958333 * u.deg * (jd_ephem - 2396758)/36525
    return Longitude(k)
    
def solar_north(t='now'):
    """Returns the position of the Solar north pole in degrees."""
    T = julian_centuries(t)
    ob1 = true_obliquity_of_ecliptic(t)
    # in degrees
    i = constants.constant('inclination solar equator')
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
    i = constants.constant('inclination solar equator')
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
