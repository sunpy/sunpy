import pytest

from astropy import units as u

from sunpy.net.hek.utils import _parse_unit


@pytest.mark.parametrize(('input_unit', 'expected_unit'), [
    ('DN/sec/pixel', u.DN/(u.pix*u.s)),
    ('ergs per cubic centimeter', u.erg/u.cm**3),
    ('HMI pixel', 0.5*u.arcsec/u.pixel),
    ('HMI pixels', 0.5*u.arcsec/u.pixel),
    ('ma/m2', u.milliampere/u.m**2),
])
def test_parse_unit(input_unit, expected_unit):
    assert _parse_unit(input_unit) == expected_unit
