from astropy import units as u

from sunpy.net.hek.utils import get_unit


def test_get_unit():
    erg_per_cm3 = u.def_unit("ergs/cm^3", u.erg/u.ml)

    u1 = get_unit('DN/sec/pixel')
    u2 = get_unit('ergs per cubic centimeter')
    u3 = get_unit('m/s/s')

    assert u1 == u.Unit('DN / (pix s)')
    assert u2 == erg_per_cm3
    assert u3 == u.Unit('m / s2')
