# -*- coding: utf-8 -*-

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import EarthLocation

from ..utilities import *


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
    # Validate against published values from the Astronomical Almanac (2013)
    assert_quantity_allclose(get_sun_B0('2013-Apr-01'), -6.54*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(get_sun_B0('2013-Dec-01'), 0.88*u.deg, atol=5e-3*u.deg)

    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(get_sun_B0('1992-Oct-13'), 5.99*u.deg, atol=5e-3*u.deg)


def test_get_sun_L0():
    # Validate against published values from the Astronomical Almanac (2013)
    assert_quantity_allclose(get_sun_L0('2013-Apr-01'), 221.44*u.deg, atol=3e-2*u.deg)
    assert_quantity_allclose(get_sun_L0('2013-Dec-01'), 237.83*u.deg, atol=3e-2*u.deg)

    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(get_sun_L0('1992-Oct-13'), 238.6317*u.deg, atol=5e-5*u.deg)


def test_get_sun_P():
    # Validate against published values from the Astronomical Almanac (2013)
    assert_quantity_allclose(get_sun_P('2013-Apr-01'), -26.15*u.deg, atol=1e-2*u.deg)
    assert_quantity_allclose(get_sun_P('2013-Dec-01'), 16.05*u.deg, atol=1e-2*u.deg)

    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(get_sun_P('1992-Oct-13'), 26.27*u.deg, atol=5e-3*u.deg)


def test_get_sunearth_distance():
    # Validate against published values from the Astronomical Almanac (2013)
    assert_quantity_allclose(get_sunearth_distance('2013-Apr-01'), 0.9992311*u.AU, atol=5e-7*u.AU)
    assert_quantity_allclose(get_sunearth_distance('2013-Dec-01'), 0.9861362*u.AU, atol=5e-7*u.AU)

    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(get_sunearth_distance('1992-Oct-13'), 0.997608*u.AU, atol=5e-7*u.AU)


def test_get_sun_orientation():
    # Not currently aware of a published value to check against, so just self-check for now

    # Check the Northern Hemisphere
    angle = get_sun_orientation(EarthLocation(lat=40*u.deg, lon=-75*u.deg), '2017-07-18 12:00')
    assert_quantity_allclose(angle, -59.4*u.deg, atol=0.1*u.deg)

    # Check the Southern Hemisphere
    angle = get_sun_orientation(EarthLocation(lat=-40*u.deg, lon=-75*u.deg), '2017-02-18 13:00')
    assert_quantity_allclose(angle, -110.8*u.deg, atol=0.1*u.deg)
