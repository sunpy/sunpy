import astropy.units as u

from sunpy.net import attr
from sunpy.net import attrs as a
from sunpy.time import parse_time


def test_time_subclass():
    class NewTime(a.Time):
        pass

    assert isinstance(NewTime("2020/01/01", "2020/01/02") & a.Time("2020/02/02", "2020/02/03"), attr.AttrAnd)


def test_attrs_time():
    times = a.Time("2020/10/01T00:00", "2020/10/01T00:00")
    assert times.start == parse_time("2020/10/01T00:00")
    assert times.end == parse_time("2020/10/01T00:00")


def test_wavelength_attr():
    # Wavelength converts the input to be either kHz, keV or Angstrom
    wave = a.Wavelength(1*u.Mm, 2*u.Mm)
    assert u.allclose(wave.min, 1.e+16*u.Angstrom)
    assert wave.min.unit == u.Angstrom
    assert u.allclose(wave.max, 2.e+16*u.Angstrom)
    assert wave.max.unit == u.Angstrom
    assert wave.unconverted_value == (1*u.Mm, 2*u.Mm)
