# -*- coding: utf-8 -*-

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from ..utilities import get_earth


def test_get_earth():
    # Validate against published values from the Astronomical Almanac (2013)
    e1 = get_earth('2013-Jan-01')
    assert_quantity_allclose(e1.lon, 0*u.deg, atol=1e-10*u.deg)  # 0 in heliographic Stonyhurst
    assert_quantity_allclose(e1.lat, -3.03*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e1.radius, 0.9832947*u.AU, atol=5e-7*u.AU)

    e2 = get_earth('2013-Sep-01')
    assert_quantity_allclose(e2.lon, 0*u.deg, atol=1e-10*u.deg)  # 0 in heliographic Stonyhurst
    assert_quantity_allclose(e2.lat, 7.19*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(e2.radius, 1.0092561*u.AU, atol=5e-7*u.AU)
