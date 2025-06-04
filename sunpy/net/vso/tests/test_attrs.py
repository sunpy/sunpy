import pytest

import astropy.units as u

from sunpy.net import _attrs as core_attrs
from sunpy.net import attr
from sunpy.net import attrs as a
from sunpy.net.vso import attrs as va
from sunpy.time import TimeRange, parse_time


def test_simpleattr_apply():
    a = attr.ValueAttr({('test', ): 1})
    dct = {}
    va._walker.apply(a, None, dct)
    assert dct['test'] == 1


def test_Time_timerange():
    t = core_attrs.Time(TimeRange('2012/1/1', '2012/1/2'))
    assert isinstance(t, core_attrs.Time)
    assert t.min == parse_time((2012, 1, 1))
    assert t.max == parse_time((2012, 1, 2))


def test_input_error():
    with pytest.raises(ValueError, match="Specify start and end or start has to be a TimeRange"):
        core_attrs.Time('2012/1/1')


@pytest.mark.remote_data
def test_simpleattr_create(client):
    a = attr.ValueAttr({('instrument', ): 'eit'})
    assert va._walker.create(a, client.api)[0].instrument == 'eit'


def test_simpleattr_and_duplicate():
    attr = core_attrs.Instrument('foo')
    pytest.raises(TypeError, lambda: attr & core_attrs.Instrument('bar'))
    attr |= a.Source('foo')
    pytest.raises(TypeError, lambda: attr & core_attrs.Instrument('bar'))
    otherattr = core_attrs.Instrument('foo') | a.Source('foo')
    pytest.raises(TypeError, lambda: attr & otherattr)
    pytest.raises(TypeError, lambda: (attr | otherattr) & core_attrs.Instrument('bar'))
    tst = core_attrs.Instrument('foo') & a.Source('foo')
    pytest.raises(TypeError, lambda: tst & tst)


def test_simpleattr_or_eq():
    attr = core_attrs.Instrument('eit')

    assert attr | attr == attr
    assert attr | core_attrs.Instrument('eit') == attr


def test_complexattr_apply():
    tst = {('test', 'foo'): 'a', ('test', 'bar'): 'b'}
    a = attr.ValueAttr(tst)
    dct = {'test': {}}
    va._walker.apply(a, None, dct)
    assert dct['test'] == {'foo': 'a', 'bar': 'b'}


@pytest.mark.remote_data
def test_complexattr_create(client):
    a = attr.ValueAttr({('time', 'start'): 'test'})
    assert va._walker.create(a, client.api)[0].time['start'] == 'test'


def test_complexattr_and_duplicate():
    attr = core_attrs.Time((2011, 1, 1), (2011, 1, 1, 1))
    pytest.raises(TypeError,
                  lambda: attr & core_attrs.Time((2011, 2, 1), (2011, 2, 1, 1)))
    attr |= a.Source('foo')
    pytest.raises(TypeError,
                  lambda: attr & core_attrs.Time((2011, 2, 1), (2011, 2, 1, 1)))


def test_complexattr_or_eq():
    attr = core_attrs.Time((2011, 1, 1), (2011, 1, 1, 1))

    assert attr | attr == attr
    assert attr | core_attrs.Time((2011, 1, 1), (2011, 1, 1, 1)) == attr


def test_attror_and():
    attr = core_attrs.Instrument('foo') | core_attrs.Instrument('bar')
    one = attr & a.Source('bar')
    other = ((core_attrs.Instrument('foo') & a.Source('bar')) |
             (core_attrs.Instrument('bar') & a.Source('bar')))
    assert one == other


def test_wave_input_quantity():
    wrong_type_message = "Wave inputs ([10, 23]) must be astropy Quantities"
    with pytest.raises(TypeError) as excinfo:
        core_attrs.Wavelength(10, 23)
    assert wrong_type_message in str(excinfo.value)
    wrong_type_message = "Wave inputs ([<Quantity 10. Angstrom>, 23]) must be astropy Quantities"
    with pytest.raises(TypeError) as excinfo:
        core_attrs.Wavelength(10 * u.AA, 23)
    assert wrong_type_message in str(excinfo.value)


def test_wave_toangstrom():
    # TODO: this test should test that inputs are in any of spectral units
    # more than just converted to Angstroms.
    frequency = [(1, 1 * u.Hz),
                 (1e3, 1 * u.kHz),
                 (1e6, 1 * u.MHz),
                 (1e9, 1 * u.GHz)]

    energy = [(1, 1 * u.eV),
              (1e3, 1 * u.keV),
              (1e6, 1 * u.MeV)]

    for factor, unit in energy:
        w = core_attrs.Wavelength((62 / factor) * unit, (62 / factor) * unit)
        assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    w = core_attrs.Wavelength(62 * u.eV, 62 * u.eV)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199
    w = core_attrs.Wavelength(62e-3 * u.keV, 62e-3 * u.keV)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    for factor, unit in frequency:
        w = core_attrs.Wavelength((1.506e16 / factor) * unit, (1.506e16 / factor) * unit)
        assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    w = core_attrs.Wavelength(1.506e16 * u.Hz, 1.506e16 * u.Hz)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199
    w = core_attrs.Wavelength(1.506e7 * u.GHz, 1.506e7 * u.GHz)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    with pytest.raises(u.UnitsError) as excinfo:
        core_attrs.Wavelength(10 * u.g, 23 * u.g)
    assert ('This unit is not convertible to any of [Unit("Angstrom"), Unit("kHz"), '
            'Unit("keV")]' in str(excinfo.value))


def test_time_xor():
    one = core_attrs.Time((2010, 1, 1), (2010, 1, 2))
    a = one ^ core_attrs.Time((2010, 1, 1, 1), (2010, 1, 1, 2))

    assert a == attr.AttrOr(
        [core_attrs.Time((2010, 1, 1), (2010, 1, 1, 1)),
         core_attrs.Time((2010, 1, 1, 2), (2010, 1, 2))])

    a ^= core_attrs.Time((2010, 1, 1, 4), (2010, 1, 1, 5))
    assert a == attr.AttrOr([
        core_attrs.Time((2010, 1, 1), (2010, 1, 1, 1)),
        core_attrs.Time((2010, 1, 1, 2), (2010, 1, 1, 4)),
        core_attrs.Time((2010, 1, 1, 5), (2010, 1, 2))
    ])


def test_wave_xor():
    one = core_attrs.Wavelength(0 * u.AA, 1000 * u.AA)
    a = one ^ core_attrs.Wavelength(200 * u.AA, 400 * u.AA)

    assert a == attr.AttrOr([core_attrs.Wavelength(0 * u.AA, 200 * u.AA),
                             core_attrs.Wavelength(400 * u.AA, 1000 * u.AA)])

    a ^= core_attrs.Wavelength(600 * u.AA, 800 * u.AA)

    assert a == attr.AttrOr(
        [core_attrs.Wavelength(0 * u.AA, 200 * u.AA), core_attrs.Wavelength(400 * u.AA, 600 * u.AA),
         core_attrs.Wavelength(800 * u.AA, 1000 * u.AA)])


def test_err_dummyattr_create():
    with pytest.raises(TypeError):
        va._walker.create(attr.DummyAttr(), None, {})


def test_err_dummyattr_apply():
    with pytest.raises(TypeError):
        va._walker.apply(attr.DummyAttr(), None, {})


def test_wave_repr():
    """Tests the __repr__ method of class vso.attrs.Wave"""
    wav = core_attrs.Wavelength(12 * u.AA, 16 * u.AA)
    moarwav = core_attrs.Wavelength(15 * u.AA, 12 * u.AA)
    assert repr(wav) == "<sunpy.net.attrs.Wavelength(12.0, 16.0, 'Angstrom')>"
    assert repr(moarwav) == "<sunpy.net.attrs.Wavelength(12.0, 15.0, 'Angstrom')>"


def test_construct_extent():
    # yes this is coverage bingo
    ext = va.Extent(10, 20, 30, 40, 'FULLDISK')
    assert ext.x == 10
    assert ext.y == 20
    assert ext.width == 30
    assert ext.length == 40
    assert ext.type == 'FULLDISK'
