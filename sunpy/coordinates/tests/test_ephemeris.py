# -*- coding: utf-8 -*-

import pytest

import astropy.units as u
from astropy.config.paths import set_temp_cache
from astropy.constants import c as speed_of_light
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from sunpy.coordinates.ephemeris import *


def test_get_body_heliographic_stonyhurst():
    # Validate against published values from the Astronomical Almanac (2013)
    e1 = get_body_heliographic_stonyhurst('earth', '2013-Jan-01')
    assert_quantity_allclose(e1.lon, 0*u.deg, atol=1e-12*u.deg)
    assert_quantity_allclose(e1.lat, -3.03*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e1.radius, 0.9832947*u.AU, atol=5e-7*u.AU)

    e2 = get_body_heliographic_stonyhurst('earth', '2013-Sep-01')
    # https://github.com/sunpy/sunpy/issues/2727
    assert_quantity_allclose((e2.lon+1*u.deg)%(360*u.deg), 1*u.deg, atol=1e-12*u.deg)
    assert_quantity_allclose(e2.lat, 7.19*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e2.radius, 1.0092561*u.AU, atol=5e-7*u.AU)


def test_get_body_heliographic_stonyhurst_light_travel_time():
    # Tests whether the apparent position of the Sun accoutns for light travel time
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


def test_get_sun_B0():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(get_sun_B0('1992-Oct-13'), 5.99*u.deg, atol=5e-3*u.deg)


def test_get_sun_B0_array_time():
    # Validate against published values from the Astronomical Almanac (2013)
    sun_B0 = get_sun_B0(Time(['2013-04-01', '2013-12-01']))
    assert_quantity_allclose(sun_B0[0], -6.54*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(sun_B0[1], 0.88*u.deg, atol=5e-3*u.deg)


def test_get_sun_L0():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(get_sun_L0('1992-Oct-13'), 238.6317*u.deg, atol=5e-5*u.deg)


def test_get_sun_L0_array_time():
    # Validate against published values from the Astronomical Almanac (2013)
    sun_L0 = get_sun_L0(Time(['2013-04-01', '2013-12-01']))
    assert_quantity_allclose(sun_L0[0], 221.44*u.deg, atol=3e-2*u.deg)
    assert_quantity_allclose(sun_L0[1], 237.83*u.deg, atol=3e-2*u.deg)


def test_get_sun_P():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(get_sun_P('1992-Oct-13'), 26.27*u.deg, atol=5e-3*u.deg)


def test_get_sun_P_array_time():
    # Validate against published values from the Astronomical Almanac (2013)
    sun_P = get_sun_P(Time(['2013-04-01', '2013-12-01']))
    assert_quantity_allclose(sun_P[0], -26.15*u.deg, atol=1e-2*u.deg)
    assert_quantity_allclose(sun_P[1], 16.05*u.deg, atol=1e-2*u.deg)


def test_get_sunearth_distance():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(get_sunearth_distance('1992-Oct-13'), 0.997608*u.AU, atol=5e-7*u.AU)


def test_get_sunearth_distance_array_time():
    # Validate against published values from the Astronomical Almanac (2013)
    sunearth_distance = get_sunearth_distance(Time(['2013-04-01', '2013-12-01']))
    assert_quantity_allclose(sunearth_distance[0], 0.9992311*u.AU, atol=5e-7*u.AU)
    assert_quantity_allclose(sunearth_distance[1], 0.9861362*u.AU, atol=5e-7*u.AU)


def test_get_sun_orientation():
    # Not currently aware of a published value to check against, so just self-check for now

    # Check the Northern Hemisphere
    angle = get_sun_orientation(EarthLocation(lat=40*u.deg, lon=-75*u.deg), '2017-07-18 12:00')
    assert_quantity_allclose(angle, -59.4*u.deg, atol=0.1*u.deg)

    # Check the Southern Hemisphere
    angle = get_sun_orientation(EarthLocation(lat=-40*u.deg, lon=-75*u.deg), '2017-02-18 13:00')
    assert_quantity_allclose(angle, -110.8*u.deg, atol=0.1*u.deg)



@pytest.mark.remote_data
def test_get_horizons_coord(tmpdir):
    # get_horizons_coord() depends on astroquery
    astroquery = pytest.importorskip("astroquery")

    with set_temp_cache(tmpdir):
        # Validate against published values from the Astronomical Almanac (2013)
        e1 = get_horizons_coord('Geocenter', '2013-Jan-01')
    assert_quantity_allclose((e1.lon + 1*u.deg) % (360*u.deg), 1*u.deg, atol=5e-6*u.deg)
    assert_quantity_allclose(e1.lat, -3.03*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e1.radius, 0.9832947*u.AU, atol=5e-7*u.AU)

    with set_temp_cache(tmpdir):
        e2 = get_horizons_coord('Geocenter', '2013-Sep-01')
    assert_quantity_allclose((e2.lon + 1*u.deg) % (360*u.deg), 1*u.deg, atol=5e-6*u.deg)
    assert_quantity_allclose(e2.lat, 7.19*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e2.radius, 1.0092561*u.AU, atol=5e-7*u.AU)


@pytest.mark.remote_data
def test_consistency_with_horizons(tmpdir):
    # get_horizons_coord() depends on astroquery
    astroquery = pytest.importorskip("astroquery")

    # Check whether the location of Earth is the same between Astropy and JPL HORIZONS
    e1 = get_earth()
    with set_temp_cache(tmpdir):
        e2 = get_horizons_coord('Geocenter')
    assert_quantity_allclose(e1.separation_3d(e2), 0*u.km, atol=25*u.km)
