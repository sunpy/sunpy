from __future__ import absolute_import

from sunpy.sun import sun
from numpy.testing import assert_array_almost_equal

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

def test_mean_anomaly():
    assert_array_almost_equal(sun.mean_anomaly("2002/12/12"), 337.538, decimal=0)
    assert_array_almost_equal(sun.mean_anomaly("2003/03/25"), 79.055, decimal=0)
    assert_array_almost_equal(sun.mean_anomaly("2005/06/05"), 150.492, decimal=0)
    assert_array_almost_equal(sun.mean_anomaly("2006/11/17"), 312.860, decimal=0)
    assert_array_almost_equal(sun.mean_anomaly("2008/07/29"), 203.933, decimal=0)
    assert_array_almost_equal(sun.mean_anomaly("2011/01/31"), 26.742, decimal=0)
