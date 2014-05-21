from __future__ import absolute_import

from sunpy.sun import sun
from numpy.testing import assert_array_almost_equal
from astropy import units as u

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

def test_mean_anomaly():
    assert_array_almost_equal(sun.mean_anomaly("2002/12/12"), 337.538 * u.deg, decimal=0)
    assert_array_almost_equal(sun.mean_anomaly("2003/03/25"), 79.055 * u.deg, decimal=0)
    assert_array_almost_equal(sun.mean_anomaly("2005/06/05"), 150.492 * u.deg, decimal=0)
    assert_array_almost_equal(sun.mean_anomaly("2006/11/17"), 312.860 * u.deg, decimal=0)
    assert_array_almost_equal(sun.mean_anomaly("2008/07/29"), 203.933 * u.deg, decimal=0)
    assert_array_almost_equal(sun.mean_anomaly("2011/01/31"), 26.742 * u.deg, decimal=0)
