# -*- coding: utf-8 -*-

import pytest

import astropy.units as u
from astropy.config.paths import set_temp_cache
from astropy.constants import c as speed_of_light
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from sunpy.coordinates.ephemeris import *
from sunpy.time import parse_time


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


@pytest.mark.remote_data
def test_get_horizons_coord():
    # get_horizons_coord() depends on astroquery
    astroquery = pytest.importorskip("astroquery")

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
    astroquery = pytest.importorskip("astroquery")

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


@pytest.mark.remote_data
def test_consistency_with_horizons():
    # get_horizons_coord() depends on astroquery
    astroquery = pytest.importorskip("astroquery")

    # Check whether the location of Earth is the same between Astropy and JPL HORIZONS
    now = parse_time('now')
    e1 = get_earth(now)
    e2 = get_horizons_coord('Geocenter', now)
    assert_quantity_allclose(e1.separation_3d(e2), 0*u.km, atol=35*u.km)
