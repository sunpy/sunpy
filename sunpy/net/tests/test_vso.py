# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

#pylint: disable=W0613

from __future__ import absolute_import

import tempfile
import datetime

import pytest
from astropy import units as u

from sunpy.time import TimeRange
from sunpy.net import vso
from sunpy.net.vso import attrs as va
from sunpy.net.vso.vso import QueryResponse
from sunpy.net import attr


@pytest.fixture
def eit(request):
    return va.Instrument('eit')


@pytest.fixture
def client(request):
    return vso.VSOClient()


@pytest.fixture
def iclient(request):
    return vso.InteractiveVSOClient()


def test_simpleattr_apply():
    a = attr.ValueAttr({('test', ): 1})
    dct = {}
    va.walker.apply(a, None, dct)
    assert dct['test'] == 1


def test_Time_timerange():
    t = va.Time(TimeRange('2012/1/1', '2012/1/2'))
    assert isinstance(t, va.Time)
    assert t.min == datetime.datetime(2012, 1, 1)
    assert t.max == datetime.datetime(2012, 1, 2)


def test_input_error():
    with pytest.raises(ValueError):
        va.Time('2012/1/1')


@pytest.mark.online
def test_simpleattr_create(client):
    a = attr.ValueAttr({('instrument', ): 'eit'})
    assert va.walker.create(a, client.api)[0].instrument == 'eit'


def test_simpleattr_and_duplicate():
    attr = va.Instrument('foo')
    pytest.raises(TypeError, lambda: attr & va.Instrument('bar'))
    attr |= va.Source('foo')
    pytest.raises(TypeError, lambda: attr & va.Instrument('bar'))
    otherattr = va.Instrument('foo') | va.Source('foo')
    pytest.raises(TypeError, lambda: attr & otherattr)
    pytest.raises(TypeError, lambda: (attr | otherattr) & va.Instrument('bar'))
    tst = va.Instrument('foo') & va.Source('foo')
    pytest.raises(TypeError, lambda: tst & tst)


def test_simpleattr_or_eq():
    attr = va.Instrument('eit')

    assert attr | attr == attr
    assert attr | va.Instrument('eit') == attr


def test_complexattr_apply():
    tst = {('test', 'foo'): 'a', ('test', 'bar'): 'b'}
    a = attr.ValueAttr(tst)
    dct = {'test': {}}
    va.walker.apply(a, None, dct)
    assert dct['test'] == {'foo': 'a', 'bar': 'b'}


@pytest.mark.online
def test_complexattr_create(client):
    a = attr.ValueAttr({('time', 'start'): 'test'})
    assert va.walker.create(a, client.api)[0].time.start == 'test'


def test_complexattr_and_duplicate():
    attr = va.Time((2011, 1, 1), (2011, 1, 1, 1))
    pytest.raises(TypeError,
                  lambda: attr & va.Time((2011, 2, 1), (2011, 2, 1, 1)))
    attr |= va.Source('foo')
    pytest.raises(TypeError,
                  lambda: attr & va.Time((2011, 2, 1), (2011, 2, 1, 1)))


def test_complexattr_or_eq():
    attr = va.Time((2011, 1, 1), (2011, 1, 1, 1))

    assert attr | attr == attr
    assert attr | va.Time((2011, 1, 1), (2011, 1, 1, 1)) == attr


def test_attror_and():
    attr = va.Instrument('foo') | va.Instrument('bar')
    one = attr & va.Source('bar')
    other = ((va.Instrument('foo') & va.Source('bar')) |
             (va.Instrument('bar') & va.Source('bar')))
    assert one == other


def test_wave_inputQuantity():
    wrong_type_mesage = "Wave inputs must be astropy Quantities"
    with pytest.raises(TypeError) as excinfo:
        va.Wavelength(10, 23)
        assert excinfo.value.message == wrong_type_mesage
    with pytest.raises(TypeError) as excinfo:
        va.Wavelength(10 * u.AA, 23)
        assert excinfo.value.message == wrong_type_mesage


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
        w = va.Wavelength((62 / factor) * unit, (62 / factor) * unit)
        assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    w = va.Wavelength(62 * u.eV, 62 * u.eV)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199
    w = va.Wavelength(62e-3 * u.keV, 62e-3 * u.keV)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    for factor, unit in frequency:
        w = va.Wavelength((1.506e16 / factor) * unit, (1.506e16 / factor) * unit)
        assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    w = va.Wavelength(1.506e16 * u.Hz, 1.506e16 * u.Hz)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199
    w = va.Wavelength(1.506e7 * u.GHz, 1.506e7 * u.GHz)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    with pytest.raises(u.UnitsError) as excinfo:
        va.Wavelength(10 * u.g, 23 * u.g)
    assert 'This unit is not convertable to any of [Unit("Angstrom"), Unit("kHz"), Unit("keV")]' in str(excinfo)


def test_time_xor():
    one = va.Time((2010, 1, 1), (2010, 1, 2))
    a = one ^ va.Time((2010, 1, 1, 1), (2010, 1, 1, 2))

    assert a == attr.AttrOr([va.Time((2010, 1, 1), (2010, 1, 1, 1)), va.Time(
        (2010, 1, 1, 2), (2010, 1, 2))])

    a ^= va.Time((2010, 1, 1, 4), (2010, 1, 1, 5))
    assert a == attr.AttrOr([va.Time((2010, 1, 1), (2010, 1, 1, 1)), va.Time(
        (2010, 1, 1, 2), (2010, 1, 1, 4)), va.Time((2010, 1, 1, 5),
                                                   (2010, 1, 2))])


def test_wave_xor():
    one = va.Wavelength(0 * u.AA, 1000 * u.AA)
    a = one ^ va.Wavelength(200 * u.AA, 400 * u.AA)

    assert a == attr.AttrOr([va.Wavelength(0 * u.AA, 200 * u.AA), va.Wavelength(400 * u.AA, 1000 * u.AA)])

    a ^= va.Wavelength(600 * u.AA, 800 * u.AA)

    assert a == attr.AttrOr(
        [va.Wavelength(0 * u.AA, 200 * u.AA), va.Wavelength(400 * u.AA, 600 * u.AA), va.Wavelength(800 * u.AA, 1000 * u.AA)])


def test_err_dummyattr_create():
    with pytest.raises(TypeError):
        va.walker.create(attr.DummyAttr(), None, {})


def test_err_dummyattr_apply():
    with pytest.raises(TypeError):
        va.walker.apply(attr.DummyAttr(), None, {})


def test_wave_repr():
    """Tests the __repr__ method of class vso.attrs.Wave"""
    wav = vso.attrs.Wavelength(12 * u.AA, 16 * u.AA)
    moarwav = vso.attrs.Wavelength(15 * u.AA, 12 * u.AA)
    assert repr(wav) == "<Wavelength(12.0, 16.0, 'Angstrom')>"
    assert repr(moarwav) == "<Wavelength(12.0, 15.0, 'Angstrom')>"

def test_str():
    qr = QueryResponse([])
    assert str(qr) == 'Start Time End Time Source Instrument Type\n---------- -------- ------ ---------- ----'

def test_repr():
    qr = QueryResponse([])
    assert "Start Time End Time  Source Instrument   Type" in repr(qr)

@pytest.mark.online
def test_path(client):
    """
    Test that '{file}' is automatically appended to the end of a custom path if
    it is not specified.
    """
    qr = client.query(va.Time('2011-06-07 06:33',
                              '2011-06-07 06:33:08'),
                      va.Instrument('aia'),
                      va.Wavelength(171*u.AA))
    tmp_dir = tempfile.mkdtemp()
    files = client.get(qr, path=tmp_dir).wait(progress=False)

    assert len(files) == 1

    # The construction of a VSO filename is bonkers complex, so there is no
    # practical way to determine what it should be in this test, so we just
    # put it here.
    assert "aia_lev1_171a_2011_06_07t06_33_02_77z_image_lev1.fits" in files[0]
