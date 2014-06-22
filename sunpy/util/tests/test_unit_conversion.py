from sunpy.util import unit_conversion as ut
from numpy.testing import assert_array_almost_equal as a
from astropy import units as u


def test_kelvin_to_keV_vice_versa():
    a(ut.keV_to_kelvin(ut.kelvin_to_keV(300 * u.K)), 300 * u.K, decimal=7)
    a(ut.keV_to_kelvin(ut.kelvin_to_keV(273.5 * u.K)), 273.5 * u.K, decimal=7)
    a(ut.keV_to_kelvin(ut.kelvin_to_keV(3023.12 * u.K)), 3023.12 * u.K,
                                                         decimal=7)
