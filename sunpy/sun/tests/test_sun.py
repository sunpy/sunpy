from astropy.coordinates import Angle
from astropy.time import Time
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from sunpy.sun import sun


def test_true_longitude():
    # Validate against a published value from the Astronomical Almanac (1992)
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.true_longitude(t), Angle('199d54m26.17s'), atol=0.1*u.arcsec)


def test_apparent_longitude():
    # Validate against a published value from the Astronomical Almanac (1992)
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.apparent_longitude(t), Angle('199d54m21.56s'), atol=0.1*u.arcsec)


def test_true_latitude():
    # Validate against a published value from the Astronomical Almanac (1992)
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.true_latitude(t), Angle('0.72s'), atol=0.05*u.arcsec)


def test_apparent_latitude():
    # Validate against a published value from the Astronomical Almanac (1992)
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.apparent_latitude(t), Angle('0.72s'), atol=0.05*u.arcsec)


def test_solar_cycle_number():
    assert_quantity_allclose(sun.solar_cycle_number("2012/11/11"), 5, atol=1e-1)
    assert_quantity_allclose(sun.solar_cycle_number("2011/2/22"), 4, atol=1e-1)
    assert_quantity_allclose(sun.solar_cycle_number("2034/1/15"), 27, atol=1e-1)


def test_solar_semidiameter_angular_size():
    assert_quantity_allclose(sun.solar_semidiameter_angular_size("2012/11/11"), 968.871294 * u.arcsec, atol=1e-3 * u.arcsec)
    assert_quantity_allclose(sun.solar_semidiameter_angular_size("2043/03/01"), 968.326347 * u.arcsec, atol=1e-3 * u.arcsec)
    assert_quantity_allclose(sun.solar_semidiameter_angular_size("2001/07/21"), 944.039007 * u.arcsec, atol=1e-3 * u.arcsec)


def test_mean_obliquity_of_ecliptic():
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.mean_obliquity_of_ecliptic(t), 84384.8*u.arcsec, atol=0.1*u.arcsec)


def test_true_rightascension():
    assert_quantity_allclose(sun.true_rightascension("2012/11/11"), 226.548*u.deg, atol=1e-3*u.deg)
    assert_quantity_allclose(sun.true_rightascension("2142/02/03"), 316.466*u.deg, atol=1e-3*u.deg)
    assert_quantity_allclose(sun.true_rightascension("2013/12/11"), 258.150*u.deg, atol=1e-3*u.deg)


def test_true_rightascension_J2000():
    # Validate against JPL HORIZONS output
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.true_rightascension(t, equinox_of_date=False),
                             Angle('13h13m53.65s'), atol=0.01*u.arcsec)


def test_true_declination():
    assert_quantity_allclose(sun.true_declination("2012/11/11"), -17.470*u.deg, atol=1e-3*u.deg)
    assert_quantity_allclose(sun.true_declination("2245/12/01"), -21.717*u.deg, atol=1e-3*u.deg)
    assert_quantity_allclose(sun.true_declination("2014/05/27"), 21.245*u.deg, atol=1e-3*u.deg)


def test_true_declination_J2000():
    # Validate against JPL HORIZONS output
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.true_declination(t, equinox_of_date=False),
                             Angle('-7d49m20.8s'), atol=0.05*u.arcsec)


def test_true_obliquity_of_ecliptic():
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.true_obliquity_of_ecliptic(t), 84384.5*u.arcsec, atol=0.1*u.arcsec)


def test_apparent_rightascension():
    # Validate against a published value from the Astronomical Almanac (1992)
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.apparent_rightascension(t), Angle('13h13m30.749s'),
                             atol=0.01*u.arcsec)


def test_apparent_rightascension_J2000():
    # Regression-only test
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.apparent_rightascension(t, equinox_of_date=False),
                             Angle('13h13m52.37s'), atol=0.01*u.arcsec)


def test_apparent_declination():
    # Validate against a published value from the Astronomical Almanac (1992)
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.apparent_declination(t), Angle('-7d47m01.74s'), atol=0.05*u.arcsec)


def test_apparent_declination_J2000():
    # Regression-only test
    t = Time('1992-10-13', scale='tdb')
    assert_quantity_allclose(sun.apparent_declination(t, equinox_of_date=False),
                             Angle('-7d49m13.09s'), atol=0.05*u.arcsec)


def test_print_params():
    # Test only for any issues with printing; accuracy is covered by other tests
    sun.print_params()
