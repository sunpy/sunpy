from __future__ import absolute_import

import pytest
from numpy.testing import assert_array_almost_equal
from itertools import product
import datetime

from astropy.coordinates import Longitude
from astropy import units as u

from sunpy.sun import sun
from sunpy.time import parse_time

def test_sunearth_distance():
    # Source for these values
    # wolframalpha.com
    # http://www.wolframalpha.com/input/?i=earth-sun+distance+on+2010%2F02%2F04
    assert_array_almost_equal(sun.sunearth_distance("2010/02/04"), 0.9858 * u.AU, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2009/04/13"), 1.003 * u.AU, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2008/06/20"), 1.016 * u.AU, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2007/08/15"), 1.013 * u.AU, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2007/10/02"), 1.001 * u.AU, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2006/12/27"), 0.9834 * u.AU, decimal=3)

def test_true_longitude():
    # source: http://www.satellite-calculations.com/Satellite/suncalc.htm
    # values are deviating a little because of lack of time parameter in
    # true_longitude function
    assert_array_almost_equal(sun.true_longitude("2002/12/23"), 270.978 * u.deg, decimal=0)
    assert_array_almost_equal(sun.true_longitude("2003/01/29"), 308.661 * u.deg, decimal=0)
    assert_array_almost_equal(sun.true_longitude("2004/05/12"), 51.617 * u.deg, decimal=0)
    assert_array_almost_equal(sun.true_longitude("2006/07/04"), 101.910 * u.deg, decimal=0)
    assert_array_almost_equal(sun.true_longitude("2007/09/16"), 172.767 * u.deg, decimal=0)
    assert_array_almost_equal(sun.true_longitude("2009/02/11"), 322.394 * u.deg, decimal=0)

def test_apparent_declination():
    assert_array_almost_equal(sun.apparent_declination("2002/12/22"), -22.964 * u.deg, decimal=0)
    assert_array_almost_equal(sun.apparent_declination("2003/1/12"), -21.743 * u.deg, decimal=0)
    assert_array_almost_equal(sun.apparent_declination("2004/02/13"), -13.478 * u.deg, decimal=0)
    assert_array_almost_equal(sun.apparent_declination("2005/12/3"), -22.152 * u.deg, decimal=0)
    assert_array_almost_equal(sun.apparent_declination("2013/02/26"), -8.547 * u.deg, decimal=0)
    assert_array_almost_equal(sun.apparent_declination("2014/05/1"), 15.141 * u.deg, decimal=0)

# From Astronomical Almanac for the Year YYYY (U.s. Nautical Almanac Office) http://asa.usno.navy.mil
# (a0 + a1 * d) * u.deg
# Last parameter is the precission expected
# where d is the interval in days from YYYY January 0, 0h TT.
almanaque_mean_anomaly = {'2009': [357.528,    0.9856003,  0],
                          '2011': [357.528,    0.9856003 , 0],
                          '2013': [356.666444, 0.98560028, 6],
                          '2014': [356.410547, 0.98560028, 6],
                          '2015': [356.154649, 0.98560028, 6]}
dates = ['{}-01-01T12:00:00', '{}-05-31T13:00:00', '{}-07-31T13:00:00', '{}-08-23T00:25:00']

@pytest.mark.parametrize("year, date", product(almanaque_mean_anomaly.keys(), dates))
def test_mean_anomaly(year, date):
    d0 = datetime.datetime(int(year) - 1, 12, 31)
    date = date.format(year)
    parameters = almanaque_mean_anomaly[year]

    number_of_days = lambda d: (parse_time(d) - d0).total_seconds()/3600./24.
    almanaque_eq = lambda d: Longitude((parameters[0] + parameters[1] * d)* u.deg)
    almanaque = almanaque_eq(number_of_days(date))
    assert_array_almost_equal(sun.mean_anomaly(date), almanaque, decimal=parameters[2])

#These values are tested from the functions after the integration of astropy.units

# Peak solar cycle dates for each solar cycle from solar cycle 1.
scdates = ['1761-06', '1769-10', '1778-05', '1787-11', '1804-12',
           '1816-03', '1829-06', '1837-02', '1847-11', '1860-07',
           '1884-01', '1893-08', '1905-10', '1917-08', '1928-06',
           '1937-05', '1947-07', '1957-11', '1989-02', '1979-11',
           '1989-10', '2000-01']

@pytest.mark.parametrize('date, cycle', zip(scdates, range(1, len(scdates)+1))) 
def test_solar_cycle_number(date, cycle):
    """
    Function to test the solar cycle number extracting the peaks from
    http://users.telenet.be/j.janssens/Engzonnecyclus.html#Overzicht
    """
    assert_array_almost_equal(sun.solar_cycle_number(date+'-01'), cycle, decimal=0)


def test_solar_semidiameter_angular_size():
    assert_array_almost_equal(sun.solar_semidiameter_angular_size("2012/11/11"), 968.383 * u.arcsec, decimal=3)
    assert_array_almost_equal(sun.solar_semidiameter_angular_size("2043/03/01"), 968.274 * u.arcsec, decimal=3)
    assert_array_almost_equal(sun.solar_semidiameter_angular_size("2001/07/21"), 943.706 * u.arcsec, decimal=3)

def test_mean_ecliptic_longitude():
    assert_array_almost_equal(sun.mean_ecliptic_longitude("2012/11/11"), 229.558 * u.deg, decimal=3)
    assert_array_almost_equal(sun.mean_ecliptic_longitude("2101/04/29"), 35.824 * u.deg, decimal=3)
    assert_array_almost_equal(sun.mean_ecliptic_longitude("2003/09/15"), 172.568 * u.deg, decimal =3)

def test_equation_of_center():
    assert_array_almost_equal(sun.equation_of_center("2012/11/11"), -1.559 * u.deg, decimal=3)
    assert_array_almost_equal(sun.equation_of_center("2014/05/27"), 1.203 * u.deg, decimal=3)
    assert_array_almost_equal(sun.equation_of_center("2134/01/01"), -0.194 * u.deg, decimal=3)

def test_true_anomaly():
    assert_array_almost_equal(sun.true_anomaly("2012/11/11"), 304.837 * u.deg, decimal=3)
    assert_array_almost_equal(sun.true_anomaly("2242/06/29"), 169.055 * u.deg, decimal=3)
    assert_array_almost_equal(sun.true_anomaly("2020/01/01"), 355.715 * u.deg, decimal=3)

def test_apparent_longitude():
    assert_array_almost_equal(sun.apparent_longitude("2012/11/11"), 227.989 * u.deg, decimal=3)
    assert_array_almost_equal(sun.apparent_longitude("2014/05/27"), 64.687 * u.deg, decimal=3)
    assert_array_almost_equal(sun.apparent_longitude("2134/02/12"), 322.053 * u.deg, decimal=3)

def test_true_obliquity_of_ecliptic():
    assert_array_almost_equal(sun.true_obliquity_of_ecliptic("2012/11/11"), 23.437 * u.deg, decimal=3)
    assert_array_almost_equal(sun.true_obliquity_of_ecliptic("2132/12/29"), 23.421 * u.deg, decimal=3)
    assert_array_almost_equal(sun.true_obliquity_of_ecliptic("2002/03/15"), 23.438 * u.deg, decimal=3)

def test_true_rightasenscion():
    assert_array_almost_equal(sun.true_rightascension("2012/11/11"), 359.318 * u.deg, decimal=3)
    assert_array_almost_equal(sun.true_rightascension("2142/02/03"), 359.328 * u.deg, decimal=3)
    assert_array_almost_equal(sun.true_rightascension("2013/12/11"), 359.102 * u.deg, decimal=3)

def test_true_declination():
    assert_array_almost_equal(sun.true_declination("2012/11/11"), -0.669 * u.deg, decimal=3)
    assert_array_almost_equal(sun.true_declination("2245/12/01"), -0.380 * u.deg, decimal=3)
    assert_array_almost_equal(sun.true_declination("2014/05/27"), 0.427 * u.deg, decimal=3)

def test_apparent_obliquity_of_ecliptic():
    assert_array_almost_equal(sun.apparent_obliquity_of_ecliptic("2012/11/11"), 23.435 * u.deg, decimal=3)
    assert_array_almost_equal(sun.apparent_obliquity_of_ecliptic("2014/05/27"), 23.438 * u.deg, decimal=3)
    assert_array_almost_equal(sun.apparent_obliquity_of_ecliptic("2412/02/26"), 23.388 * u.deg, decimal=3)

def test_apparent_rightascension():
    assert_array_almost_equal(sun.apparent_rightascension("2012/11/11"), 15.035 * u.hourangle, decimal=3)
    assert_array_almost_equal(sun.apparent_rightascension("2013/12/13"), 17.282 * u.hourangle, decimal=3)
    assert_array_almost_equal(sun.apparent_rightascension("2512/04/09"), 1.134 * u.hourangle, decimal=3)

@pytest.mark.thisone
def test_longitude_ascending_node():
    # Example 29.a from Astronomical Algorithms 2nd Ed. - Jean Meeus 2005
    dateinput = "1992-10-13T00:00:00"
    kvalue = 75.6597 * u.deg
    assert_array_almost_equal(sun.longitude_ascending_node(dateinput), kvalue , decimal=4)

def test_solar_north():
    assert_array_almost_equal(sun.solar_north("2012/11/11"), 15.149 * u.deg, decimal=3)
    assert_array_almost_equal(sun.solar_north("2019/10/10"), -1.693 * u.deg, decimal=3)
    assert_array_almost_equal(sun.solar_north("2542/02/20"), 41.351 * u.deg, decimal=3)

