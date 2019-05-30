import pytest

from astropy.coordinates import Angle
from astropy.time import Time
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from sunpy.sun import sun


@pytest.fixture
def t1():
    return Time('1992-10-13', scale='tdb')


@pytest.fixture
def t2():
    return Time('2013-06-17', scale='tt')


def test_true_longitude(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C16)
    assert_quantity_allclose(sun.true_longitude(t1), Angle('199d54m26.17s'), atol=0.005*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C12)
    assert_quantity_allclose(sun.true_longitude(t2), Angle('85d58m57.46s'), atol=0.005*u.arcsec)


def test_apparent_longitude(t1, t2):
    aberration_1au = Angle('20.496s')
    # Validate against a published value from the Astronomical Almanac (1992, C2+C16+B30)
    assert_quantity_allclose(sun.apparent_longitude(t1),
                             Angle('199d54m26.17s') + Angle('15.908s') - aberration_1au / 0.9976085,
                             atol=0.005*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C2+C12+B61)
    assert_quantity_allclose(sun.apparent_longitude(t2),
                             Angle('85d58m57.46') + Angle('12.0489s') - aberration_1au / 1.0159149,
                             atol=0.005*u.arcsec)


def test_true_latitude(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C16)
    # Disagreement with published precision may be due to changes in definitions after 1992
    assert_quantity_allclose(sun.true_latitude(t1), Angle('0.72s'), atol=0.05*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C12)
    assert_quantity_allclose(sun.true_latitude(t2), Angle('-0.42s'), atol=0.005*u.arcsec)


def test_apparent_latitude(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C2+C16)
    # Disagreement with published precision may be due to changes in definitions after 1992
    assert_quantity_allclose(sun.apparent_latitude(t1), Angle('0.72s'), atol=0.05*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C2+C12)
    assert_quantity_allclose(sun.apparent_latitude(t2), Angle('-0.42s'), atol=0.005*u.arcsec)


def test_solar_semidiameter_angular_size():
    # Regression-only test
    assert_quantity_allclose(sun.solar_semidiameter_angular_size("2012/11/11"), 968.871294 * u.arcsec, atol=1e-3 * u.arcsec)
    assert_quantity_allclose(sun.solar_semidiameter_angular_size("2043/03/01"), 968.326347 * u.arcsec, atol=1e-3 * u.arcsec)
    assert_quantity_allclose(sun.solar_semidiameter_angular_size("2001/07/21"), 944.039007 * u.arcsec, atol=1e-3 * u.arcsec)


def test_mean_obliquity_of_ecliptic(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C1)
    # Note that the publication date pre-dates the IAU 2006 definition of obliquity
    assert_quantity_allclose(sun.mean_obliquity_of_ecliptic(t1), 84384.82*u.arcsec,
                             atol=0.05*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C1)
    assert_quantity_allclose(sun.mean_obliquity_of_ecliptic(t2), 84375.098*u.arcsec,
                             atol=0.005*u.arcsec)


def test_true_rightascension():
    # Regression-only test
    assert_quantity_allclose(sun.true_rightascension("2012/11/11"), 226.548*u.deg, atol=1e-3*u.deg)
    assert_quantity_allclose(sun.true_rightascension("2142/02/03"), 316.466*u.deg, atol=1e-3*u.deg)
    assert_quantity_allclose(sun.true_rightascension("2013/12/11"), 258.150*u.deg, atol=1e-3*u.deg)


def test_true_rightascension_J2000(t1, t2):
    # Validate against JPL HORIZONS output
    # https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE=OBSERVER&COMMAND=10&CENTER=500&QUANTITIES=1&ANG_FORMAT=HMS&EXTRA_PREC=YES&TIME_TYPE=TT&TLIST=2448908.5%0A2456460.5
    assert_quantity_allclose(sun.true_rightascension(t1, equinox_of_date=False),
                             Angle('13h13m53.65s'), atol=0.005*u.arcsec)
    assert_quantity_allclose(sun.true_rightascension(t2, equinox_of_date=False),
                             Angle('05h41m40.32s'), atol=0.005*u.arcsec)


def test_true_declination():
    # Regression-only test
    assert_quantity_allclose(sun.true_declination("2012/11/11"), -17.470*u.deg, atol=1e-3*u.deg)
    assert_quantity_allclose(sun.true_declination("2245/12/01"), -21.717*u.deg, atol=1e-3*u.deg)
    assert_quantity_allclose(sun.true_declination("2014/05/27"), 21.245*u.deg, atol=1e-3*u.deg)


def test_true_declination_J2000(t1, t2):
    # Validate against JPL HORIZONS output
    # https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE=OBSERVER&COMMAND=10&CENTER=500&QUANTITIES=1&ANG_FORMAT=HMS&EXTRA_PREC=YES&TIME_TYPE=TT&TLIST=2448908.5%0A2456460.5
    assert_quantity_allclose(sun.true_declination(t1, equinox_of_date=False),
                             Angle('-7d49m20.84s'), atol=0.005*u.arcsec)
    assert_quantity_allclose(sun.true_declination(t2, equinox_of_date=False),
                             Angle('23d22m13.97s'), atol=0.005*u.arcsec)


def test_true_obliquity_of_ecliptic(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C1+B30)
    # Note that the publication date pre-dates the IAU 2006 definition of obliquity
    assert_quantity_allclose(sun.true_obliquity_of_ecliptic(t1), (84384.82 - 0.308)*u.arcsec,
                             atol=0.05*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C1+B61)
    assert_quantity_allclose(sun.true_obliquity_of_ecliptic(t2), (84375.098 - 7.0153)*u.arcsec,
                             atol=0.005*u.arcsec)


def test_apparent_rightascension(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C16)
    assert_quantity_allclose(sun.apparent_rightascension(t1), Angle('13h13m30.75s'),
                             atol=0.005*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C12)
    assert_quantity_allclose(sun.apparent_rightascension(t2), Angle('5h42m28.88s'),
                             atol=0.01*u.arcsec)  # slight disagreement with published precision


def test_apparent_rightascension_J2000(t1):
    # Regression-only test
    assert_quantity_allclose(sun.apparent_rightascension(t1, equinox_of_date=False),
                             Angle('13h13m52.37s'), atol=0.005*u.arcsec)


def test_apparent_declination(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C16)
    assert_quantity_allclose(sun.apparent_declination(t1), Angle('-7d47m01.7s'), atol=0.05*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C12)
    assert_quantity_allclose(sun.apparent_declination(t2), Angle('23d22m27.8s'), atol=0.05*u.arcsec)


def test_apparent_declination_J2000(t1):
    # Regression-only test
    assert_quantity_allclose(sun.apparent_declination(t1, equinox_of_date=False),
                             Angle('-7d49m13.1s'), atol=0.05*u.arcsec)


def test_position(t1, t2):
    pos1 = sun.position(t1)
    ra1 = sun.apparent_rightascension(t1)
    dec1 = sun.apparent_declination(t1)
    assert_quantity_allclose(pos1, (ra1, dec1))

    pos2 = sun.position(t2, equinox_of_date=False)
    ra2 = sun.apparent_rightascension(t2, equinox_of_date=False)
    dec2 = sun.apparent_declination(t2, equinox_of_date=False)
    assert_quantity_allclose(pos2, (ra2, dec2))


def test_print_params():
    # Test only for any issues with printing; accuracy is covered by other tests
    sun.print_params()
