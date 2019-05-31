# -*- coding: utf-8 -*-

import pytest

import astropy.units as u
from astropy.config.paths import set_temp_cache
from astropy.constants import c as speed_of_light
from astropy.coordinates import SkyCoord
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
