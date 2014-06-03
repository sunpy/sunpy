from sunpy.util import unit_conversion as ut
from numpy.testing import assert_array_almost_equal as a
from astropy import units as u


def test_kelvin_to_keV_vice_versa():
    a(ut.keV_to_kelvin(ut.kelvin_to_keV(300 * u.K)), 300 * u.K, decimal=7)
    a(ut.keV_to_kelvin(ut.kelvin_to_keV(273.5 * u.K)), 273.5 * u.K, decimal=7)
    a(ut.keV_to_kelvin(ut.kelvin_to_keV(3023.12 * u.K)), 3023.12 * u.K,
                                                         decimal=7)


def test_degrees_to_hours():
    assert ut.degrees_to_hours(180 * u.deg) == [12.0 * u.hourangle,
                                               0.0 * u.arcmin,
                                               0.0 * u.arcsec]
    assert ut.degrees_to_hours(90.5 * u.deg) == [6.0 * u.hourangle,
                                                1.0 * u.arcmin,
                                                59.999999999999574 * u.arcsec]
    assert ut.degrees_to_hours(11 * u.deg) == [0.0 * u.hourangle,
                                              44.0 * u.arcmin,
                                              0.0 * u.arcsec]


def test_degrees_to_arc():
    assert ut.degrees_to_arc(180 * u.deg) == [180.0 * u.deg,
                                             0.0 * u.arcmin,
                                             0.0 * u.arcsec]
    assert ut.degrees_to_arc(78.12 * u.deg) == [78.0 * u.deg, 7.0 * u.arcmin,
                                               12.000000000016371 * u.arcsec]
    assert ut.degrees_to_arc(284.982 * u.deg) == [284.0 * u.deg,
                                                 58.0 * u.arcmin,
                                                 55.20000000009986 * u.arcsec]

