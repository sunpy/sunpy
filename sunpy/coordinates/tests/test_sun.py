import pytest

from astropy.coordinates import Angle, EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from sunpy.coordinates import sun


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


def test_angular_radius():
    # Regression-only test
    assert_quantity_allclose(sun.angular_radius("2012/11/11"), 968.871*u.arcsec, atol=1e-3*u.arcsec)
    assert_quantity_allclose(sun.angular_radius("2043/03/01"), 968.326*u.arcsec, atol=1e-3*u.arcsec)
    assert_quantity_allclose(sun.angular_radius("2001/07/21"), 944.039*u.arcsec, atol=1e-3*u.arcsec)


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


def test_sky_position(t1, t2):
    pos1 = sun.sky_position(t1)
    ra1 = sun.apparent_rightascension(t1)
    dec1 = sun.apparent_declination(t1)
    assert_quantity_allclose(pos1, (ra1, dec1))

    pos2 = sun.sky_position(t2, equinox_of_date=False)
    ra2 = sun.apparent_rightascension(t2, equinox_of_date=False)
    dec2 = sun.apparent_declination(t2, equinox_of_date=False)
    assert_quantity_allclose(pos2, (ra2, dec2))


def test_print_params():
    # Test only for any issues with printing; accuracy is covered by other tests
    sun.print_params()


def test_B0():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(sun.B0('1992-Oct-13'), 5.99*u.deg, atol=5e-3*u.deg)


def test_B0_array_time():
    # Validate against published values from the Astronomical Almanac (2013)
    sun_B0 = sun.B0(Time(['2013-04-01', '2013-12-01']))
    assert_quantity_allclose(sun_B0[0], -6.54*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(sun_B0[1], 0.88*u.deg, atol=5e-3*u.deg)


def test_L0():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(sun.L0('1992-Oct-13'), 238.6317*u.deg, atol=5e-5*u.deg)


def test_L0_array_time():
    # Validate against published values from the Astronomical Almanac (2013)
    sun_L0 = sun.L0(Time(['2013-04-01', '2013-12-01']))
    assert_quantity_allclose(sun_L0[0], 221.44*u.deg, atol=3e-2*u.deg)
    assert_quantity_allclose(sun_L0[1], 237.83*u.deg, atol=3e-2*u.deg)


def test_P():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(sun.P('1992-Oct-13'), 26.27*u.deg, atol=5e-3*u.deg)


def test_P_array_time():
    # Validate against published values from the Astronomical Almanac (2013)
    sun_P = sun.P(Time(['2013-04-01', '2013-12-01']))
    assert_quantity_allclose(sun_P[0], -26.15*u.deg, atol=1e-2*u.deg)
    assert_quantity_allclose(sun_P[1], 16.05*u.deg, atol=1e-2*u.deg)


def test_earth_distance():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(sun.earth_distance('1992-Oct-13'), 0.997608*u.AU, atol=5e-7*u.AU)


def test_earth_distance_array_time():
    # Validate against published values from the Astronomical Almanac (2013)
    sunearth_distance = sun.earth_distance(Time(['2013-04-01', '2013-12-01']))
    assert_quantity_allclose(sunearth_distance[0], 0.9992311*u.AU, atol=5e-7*u.AU)
    assert_quantity_allclose(sunearth_distance[1], 0.9861362*u.AU, atol=5e-7*u.AU)


def test_orientation():
    # Not currently aware of a published value to check against, so just self-check for now

    # Check the Northern Hemisphere
    angle = sun.orientation(EarthLocation(lat=40*u.deg, lon=-75*u.deg), '2017-07-18 12:00')
    assert_quantity_allclose(angle, -59.4*u.deg, atol=0.1*u.deg)

    # Check the Southern Hemisphere
    angle = sun.orientation(EarthLocation(lat=-40*u.deg, lon=-75*u.deg), '2017-02-18 13:00')
    assert_quantity_allclose(angle, -110.8*u.deg, atol=0.1*u.deg)
