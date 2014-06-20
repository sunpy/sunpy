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
    assert_array_almost_equal(sun.sunearth_distance("2010/02/04"), 0.9858, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2009/04/13"), 1.003, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2008/06/20"), 1.016, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2007/08/15"), 1.013, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2007/10/02"), 1.001, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2006/12/27"), 0.9834, decimal=3)

def test_true_longitude():
    # source: http://www.satellite-calculations.com/Satellite/suncalc.htm
    # values are deviating a little because of lack of time parameter in
    # true_longitude function
    assert_array_almost_equal(sun.true_longitude("2002/12/23"), 270.978, decimal=0)
    assert_array_almost_equal(sun.true_longitude("2003/01/29"), 308.661, decimal=0)
    assert_array_almost_equal(sun.true_longitude("2004/05/12"), 51.617, decimal=0)
    assert_array_almost_equal(sun.true_longitude("2006/07/04"), 101.910, decimal=0)
    assert_array_almost_equal(sun.true_longitude("2007/09/16"), 172.767, decimal=0)
    assert_array_almost_equal(sun.true_longitude("2009/02/11"), 322.394, decimal=0)

def test_apparent_declination():
    assert_array_almost_equal(sun.apparent_declination("2002/12/22"), -22.964, decimal=0)
    assert_array_almost_equal(sun.apparent_declination("2003/1/12"), -21.743, decimal=0)
    assert_array_almost_equal(sun.apparent_declination("2004/02/13"), -13.478, decimal=0)
    assert_array_almost_equal(sun.apparent_declination("2005/12/3"), -22.152, decimal=0)
    assert_array_almost_equal(sun.apparent_declination("2013/02/26"), -8.547, decimal=0)
    assert_array_almost_equal(sun.apparent_declination("2014/05/1"), 15.141, decimal=0)

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

