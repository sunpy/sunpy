from __future__ import absolute_import

import astropy.units as u

from sunpy.sun import sun
from sunpy.tests.helpers import assert_quantity_allclose


def test_sunearth_distance():
    # Source for these values
    # wolframalpha.com
    # http://www.wolframalpha.com/input/?i=earth-sun+distance+on+2010%2F02%2F04
    assert_quantity_allclose(sun.sunearth_distance("2010/02/04"), 0.9858 * u.AU, atol=1e-3 * u.AU)
    assert_quantity_allclose(sun.sunearth_distance("2009/04/13"), 1.003 * u.AU, atol=1e-3 * u.AU)
    assert_quantity_allclose(sun.sunearth_distance("2008/06/20"), 1.016 * u.AU, atol=1e-3 * u.AU)
    assert_quantity_allclose(sun.sunearth_distance("2007/08/15"), 1.013 * u.AU, atol=1e-3 * u.AU)
    assert_quantity_allclose(sun.sunearth_distance("2007/10/02"), 1.001 * u.AU, atol=1e-3 * u.AU)
    assert_quantity_allclose(sun.sunearth_distance("2006/12/27"), 0.9834 * u.AU, atol=1e-3 * u.AU)


def test_true_longitude():
    # source: http://www.satellite-calculations.com/Satellite/suncalc.htm
    # values are deviating a little because of lack of time parameter in
    # true_longitude function
    assert_quantity_allclose(sun.true_longitude("2002/12/23"), 270.978 * u.deg, atol=1.1  * u.deg)
    assert_quantity_allclose(sun.true_longitude("2003/01/29"), 308.661 * u.deg, atol=1.1 * u.deg)
    assert_quantity_allclose(sun.true_longitude("2004/05/12"), 51.617 * u.deg, atol=1.1 * u.deg)
    assert_quantity_allclose(sun.true_longitude("2006/07/04"), 101.910 * u.deg, atol=1.1 * u.deg)
    assert_quantity_allclose(sun.true_longitude("2007/09/16"), 172.767 * u.deg, atol=1.1 * u.deg)
    assert_quantity_allclose(sun.true_longitude("2009/02/11"), 322.394 * u.deg, atol=1.1 * u.deg)


def test_apparent_declination():
    assert_quantity_allclose(sun.apparent_declination("2002/12/22"), -22.964 * u.deg, atol=1 * u.deg)
    assert_quantity_allclose(sun.apparent_declination("2003/1/12"), -21.743 * u.deg, atol=1 * u.deg)
    assert_quantity_allclose(sun.apparent_declination("2004/02/13"), -13.478 * u.deg, atol=1 * u.deg)
    assert_quantity_allclose(sun.apparent_declination("2005/12/3"), -22.152 * u.deg, atol=1 * u.deg)
    assert_quantity_allclose(sun.apparent_declination("2013/02/26"), -8.547 * u.deg, atol=1 * u.deg)
    assert_quantity_allclose(sun.apparent_declination("2014/05/1"), 15.141 * u.deg, atol=1 * u.deg)


def test_mean_anomaly():
    assert_quantity_allclose(sun.mean_anomaly("2002/12/12"), 337.538 * u.deg, atol=1 * u.deg)
    assert_quantity_allclose(sun.mean_anomaly("2003/03/25"), 79.055 * u.deg, atol=1 * u.deg)
    assert_quantity_allclose(sun.mean_anomaly("2005/06/05"), 150.492 * u.deg, atol=1 * u.deg)
    assert_quantity_allclose(sun.mean_anomaly("2006/11/17"), 312.860 * u.deg, atol=1 * u.deg)
    assert_quantity_allclose(sun.mean_anomaly("2008/07/29"), 203.933 * u.deg, atol=1 * u.deg)
    assert_quantity_allclose(sun.mean_anomaly("2011/01/31"), 26.742 * u.deg, atol=1 * u.deg)

#These values are tested from the functions after the integration of astropy.units


def test_solar_cycle_number():
    assert_quantity_allclose(sun.solar_cycle_number("2012/11/11"), 5, atol=1e-1)
    #76
    assert_quantity_allclose(sun.solar_cycle_number("2011/2/22"), 4, atol=1e-1)
    #23
    assert_quantity_allclose(sun.solar_cycle_number("2034/1/15"), 27, atol=1e-1)


def test_solar_semidiameter_angular_size():
    assert_quantity_allclose(sun.solar_semidiameter_angular_size("2012/11/11"), 968.383 * u.arcsec, atol=1e-3 * u.arcsec)
    assert_quantity_allclose(sun.solar_semidiameter_angular_size("2043/03/01"), 968.274 * u.arcsec, atol=1e-3 * u.arcsec)
    assert_quantity_allclose(sun.solar_semidiameter_angular_size("2001/07/21"), 943.706 * u.arcsec, atol=1e-3 * u.arcsec)


def test_mean_ecliptic_longitude():
    assert_quantity_allclose(sun.mean_ecliptic_longitude("2012/11/11"), 229.558 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.mean_ecliptic_longitude("2101/04/29"), 35.824 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.mean_ecliptic_longitude("2003/09/15"), 172.568 * u.deg, atol=3 * u.deg)


def test_equation_of_center():
    assert_quantity_allclose(sun.equation_of_center("2012/11/11"), -1.559 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.equation_of_center("2014/05/27"), 1.203 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.equation_of_center("2134/01/01"), -0.194 * u.deg, atol=1e-3 * u.deg)


def test_true_anomaly():
    assert_quantity_allclose(sun.true_anomaly("2012/11/11"), 304.837 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.true_anomaly("2242/06/29"), 169.055 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.true_anomaly("2020/01/01"), 355.715 * u.deg, atol=1e-3 * u.deg)


def test_apparent_longitude():
    assert_quantity_allclose(sun.apparent_longitude("2012/11/11"), 227.989 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.apparent_longitude("2014/05/27"), 64.687 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.apparent_longitude("2134/02/12"), 322.053 * u.deg, atol=1e-3 * u.deg)


def test_true_obliquity_of_ecliptic():
    assert_quantity_allclose(sun.true_obliquity_of_ecliptic("2012/11/11"), 23.437 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.true_obliquity_of_ecliptic("2132/12/29"), 23.421 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.true_obliquity_of_ecliptic("2002/03/15"), 23.438 * u.deg, atol=1e-3 * u.deg)


def test_true_rightasenscion():
    assert_quantity_allclose(sun.true_rightascension("2012/11/11"), 359.318 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.true_rightascension("2142/02/03"), 359.328 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.true_rightascension("2013/12/11"), 359.102 * u.deg, atol=1e-3 * u.deg)


def test_true_declination():
    assert_quantity_allclose(sun.true_declination("2012/11/11"), -0.669 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.true_declination("2245/12/01"), -0.380 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.true_declination("2014/05/27"), 0.427 * u.deg, atol=1e-3 * u.deg)


def test_apparent_obliquity_of_ecliptic():
    assert_quantity_allclose(sun.apparent_obliquity_of_ecliptic("2012/11/11"), 23.435 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.apparent_obliquity_of_ecliptic("2014/05/27"), 23.438 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.apparent_obliquity_of_ecliptic("2412/02/26"), 23.388 * u.deg, atol=1e-3 * u.deg)


def test_apparent_rightascension():
    assert_quantity_allclose(sun.apparent_rightascension("2012/11/11"), 15.035 * u.hourangle, atol=1e-3 * u.hourangle)
    assert_quantity_allclose(sun.apparent_rightascension("2013/12/13"), 17.282 * u.hourangle, atol=1e-3 * u.hourangle)
    assert_quantity_allclose(sun.apparent_rightascension("2512/04/09"), 1.134 * u.hourangle, atol=1e-3 * u.hourangle)


def test_solar_north():
    assert_quantity_allclose(sun.solar_north("2012/11/11"), 15.149 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.solar_north("2019/10/10"), -1.693 * u.deg, atol=1e-3 * u.deg)
    assert_quantity_allclose(sun.solar_north("2542/02/20"), 41.351 * u.deg, atol=1e-3 * u.deg)

