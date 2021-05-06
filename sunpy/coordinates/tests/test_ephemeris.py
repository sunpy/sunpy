
import pytest
from hypothesis import HealthCheck, given, settings

import astropy.units as u
from astropy.constants import c as speed_of_light
from astropy.coordinates import SkyCoord, solar_system_ephemeris
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from sunpy.coordinates.ephemeris import get_body_heliographic_stonyhurst, get_earth, get_horizons_coord
from .strategies import times


def test_get_body_heliographic_stonyhurst():
    # Validate against published values from the Astronomical Almanac (2013)
    e1 = get_body_heliographic_stonyhurst('earth', '2013-Jan-01')
    assert_quantity_allclose(e1.lon, 0*u.deg, atol=1e-12*u.deg)
    assert_quantity_allclose(e1.lat, -3.03*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e1.radius, 0.9832947*u.AU, atol=5e-7*u.AU)

    e2 = get_body_heliographic_stonyhurst('earth', '2013-Sep-01')
    assert_quantity_allclose(e2.lon, 0*u.deg, atol=1e-12*u.deg)
    assert_quantity_allclose(e2.lat, 7.19*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e2.radius, 1.0092561*u.AU, atol=5e-7*u.AU)


def test_get_body_heliographic_stonyhurst_light_travel_time():
    # Tests whether the apparent position of the Sun accounts for light travel time
    t = Time('2012-06-05 22:34:48.350')  # An arbitrary test time

    # Use the implemented correction for light travel time
    implementation = get_body_heliographic_stonyhurst('sun', t, observer=get_earth(t))
    implementation_icrs = SkyCoord(implementation).icrs.cartesian

    # Use a manual correction for light travel time
    light_travel_time = get_earth(t).radius / speed_of_light
    manual = get_body_heliographic_stonyhurst('sun', t - light_travel_time)
    manual_icrs = SkyCoord(manual).icrs.cartesian

    difference = (implementation_icrs - manual_icrs).norm()
    assert_quantity_allclose(difference, 0*u.m, atol=1*u.m)


def test_get_body_heliographic_stonyhurst_light_travel_time_array():
    # Tests whether requesting an array of locations returns the same answers as individually
    t1 = Time('2001-02-03 04:05:06')
    t2 = Time('2011-12-13 14:15:16')

    venus1 = get_body_heliographic_stonyhurst('venus', t1, observer=get_earth(t1))
    venus2 = get_body_heliographic_stonyhurst('venus', t2, observer=get_earth(t2))
    both = get_body_heliographic_stonyhurst('venus', [t1, t2], observer=get_earth([t1, t2]))

    assert_quantity_allclose(venus1.lon, both[0].lon)
    assert_quantity_allclose(venus1.lat, both[0].lat)
    assert_quantity_allclose(venus1.radius, both[0].radius)
    assert_quantity_allclose(venus2.lon, both[1].lon)
    assert_quantity_allclose(venus2.lat, both[1].lat)
    assert_quantity_allclose(venus2.radius, both[1].radius)


def test_get_earth():
    # Validate against published values from the Astronomical Almanac (2013)
    e1 = get_earth('2013-Jan-01')
    assert e1.lon == 0*u.deg
    assert_quantity_allclose(e1.lat, -3.03*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e1.radius, 0.9832947*u.AU, atol=5e-7*u.AU)

    e2 = get_earth('2013-Sep-01')
    assert e2.lon == 0*u.deg
    assert_quantity_allclose(e2.lat, 7.19*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e2.radius, 1.0092561*u.AU, atol=5e-7*u.AU)


@pytest.mark.remote_data
def test_get_horizons_coord():
    # get_horizons_coord() depends on astroquery
    pytest.importorskip("astroquery")

    # Validate against published values from the Astronomical Almanac (2013)
    e1 = get_horizons_coord('Geocenter', '2013-Jan-01')
    assert_quantity_allclose(e1.lon, 0*u.deg, atol=5e-6*u.deg)
    assert_quantity_allclose(e1.lat, -3.03*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e1.radius, 0.9832947*u.AU, atol=5e-7*u.AU)

    e2 = get_horizons_coord('Geocenter', '2013-Sep-01')
    assert_quantity_allclose(e1.lon, 0*u.deg, atol=5e-6*u.deg)
    assert_quantity_allclose(e2.lat, 7.19*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e2.radius, 1.0092561*u.AU, atol=5e-7*u.AU)


@pytest.mark.remote_data
def test_get_horizons_coord_array_time():
    # get_horizons_coord() depends on astroquery
    pytest.importorskip("astroquery")

    # Validate against published values from the Astronomical Almanac (2013, C8-C13)
    array_time = Time(['2013-05-01', '2013-06-01', '2013-04-01', '2013-03-01'])
    e = get_horizons_coord('Geocenter', array_time)

    assert_quantity_allclose(e[0].lon, 0*u.deg, atol=5e-6*u.deg)
    assert_quantity_allclose(e[0].lat, -4.17*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e[0].radius, 1.0075271*u.AU, atol=5e-7*u.AU)

    assert_quantity_allclose(e[1].lon, 0*u.deg, atol=5e-6*u.deg)
    assert_quantity_allclose(e[1].lat, -0.66*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e[1].radius, 1.0140013*u.AU, atol=5e-7*u.AU)

    assert_quantity_allclose(e[2].lon, 0*u.deg, atol=5e-6*u.deg)
    assert_quantity_allclose(e[2].lat, -6.54*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e[2].radius, 0.9992311*u.AU, atol=5e-7*u.AU)

    assert_quantity_allclose(e[3].lon, 0*u.deg, atol=5e-6*u.deg)
    assert_quantity_allclose(e[3].lat, -7.22*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e[3].radius, 0.9908173*u.AU, atol=5e-7*u.AU)


@pytest.fixture
def use_DE440s():
    # This class is for test functions that need the Astropy ephemeris to be set to DE432s
    pytest.importorskip("astroquery")

    old_ephemeris = solar_system_ephemeris.get()
    try:
        solar_system_ephemeris.set('de440s')
    except ValueError:
        pytest.skip("The installed version of Astropy cannot set the ephemeris to DE440s")

    yield

    solar_system_ephemeris.set(old_ephemeris)


@pytest.mark.remote_data
@given(obstime=times())
@settings(deadline=5000, max_examples=10, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_consistency_with_horizons(use_DE440s, obstime):
    # Check that the high-accuracy Astropy ephemeris has been set
    assert solar_system_ephemeris.get() == 'de440s'

    # Check whether the location of Earth is the same between Astropy and JPL HORIZONS
    e1 = get_earth(obstime)
    e2 = get_horizons_coord('Geocenter', obstime)
    assert_quantity_allclose(e2.separation_3d(e1), 0*u.km, atol=50*u.m)

    # Check whether the location of Mars is the same between Astropy and JPL HORIZONS
    e1 = get_body_heliographic_stonyhurst('mars', obstime)
    e2 = get_horizons_coord('Mars barycenter', obstime)
    assert_quantity_allclose(e2.separation_3d(e1), 0*u.km, atol=500*u.m)
