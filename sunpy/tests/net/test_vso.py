import pytest

from sunpy.net import vso
from sunpy.net import attr

from sunpy.util.util import energy, frequency, wavelength

def pytest_funcarg__eit(request):
    return vso.Instrument('eit')


def pytest_funcarg__client(request):
    return vso.VSOClient()


def pytest_funcarg__iclient(request):
    return vso.InteractiveVSOClient()


def test_simpleattr_apply():
    attr = vso.ValueAttr({('test', ): 1})
    dct = {}
    vso.walker.apply(attr, None, dct)
    assert dct['test'] == 1


def test_simpleattr_create(client):
    attr = vso.ValueAttr({('instrument', ): 'eit'})
    assert vso.walker.create(attr, client.api)[0].instrument == 'eit'


def test_simpleattr_and_duplicate():
    attr = vso.Instrument('foo')
    pytest.raises(TypeError, lambda: attr & vso.Instrument('bar'))
    attr |= vso.Source('foo')
    pytest.raises(TypeError, lambda: attr & vso.Instrument('bar'))
    otherattr = vso.Instrument('foo') | vso.Source('foo')
    pytest.raises(TypeError, lambda: attr & otherattr)
    pytest.raises(TypeError, lambda: (attr | otherattr) & vso.Instrument('bar'))
    tst = vso.Instrument('foo') & vso.Source('foo')
    pytest.raises(TypeError, lambda: tst & tst)


def test_simpleattr_or_eq():
    attr = vso.Instrument('eit')
    
    assert attr | attr == attr
    assert attr | vso.Instrument('eit') == attr


def test_complexattr_apply():
    tst = {('test', 'foo'): 'a', ('test', 'bar'): 'b'}
    attr = vso.ValueAttr(tst)
    dct = {'test': {}}
    vso.walker.apply(attr, None, dct)
    assert dct['test'] == {'foo': 'a', 'bar': 'b'}


def test_complexattr_create(client):
    attr = vso.ValueAttr({('time', 'start'): 'test'})
    assert vso.walker.create(attr, client.api)[0].time.start == 'test'


def test_complexattr_and_duplicate():
    attr = vso.Time.dt((2011, 1, 1), (2011, 1, 1, 1))
    pytest.raises(
        TypeError,
        lambda: attr & vso.Time.dt((2011, 2, 1), (2011, 2, 1, 1))
    )
    attr |= vso.Source('foo')
    pytest.raises(
        TypeError,
        lambda: attr & vso.Time.dt((2011, 2, 1), (2011, 2, 1, 1))
    )


def test_complexattr_or_eq():
    attr = vso.Time.dt((2011, 1, 1), (2011, 1, 1, 1))
    
    assert attr | attr == attr
    assert attr | vso.Time.dt((2011, 1, 1), (2011, 1, 1, 1)) == attr


def test_attror_and():
    attr = vso.Instrument('foo') | vso.Instrument('bar')
    one = attr & vso.Source('bar')
    other = (
        (vso.Instrument('foo') & vso.Source('bar')) | 
        (vso.Instrument('bar') & vso.Source('bar'))
    )
    assert one == other


def test_wave_toangstrom():
    for name, factor in energy:
        w = vso.Wave(62 / factor, 62 / factor, name)
        assert int(w.min) == 199
    
    w = vso.Wave(62, 62, 'eV')
    assert int(w.min) == 199
    w = vso.Wave(62e-3, 62e-3, 'keV')
    assert int(w.min) == 199

    for name, factor in frequency:
        w = vso.Wave(1.506e16 / factor, 1.506e16 / factor, name)
        assert int(w.min) == 199
    
    w = vso.Wave(1.506e16, 1.506e16, 'Hz')
    assert int(w.min) == 199
    w = vso.Wave(1.506e7, 1.506e7, 'GHz')
    assert int(w.min) == 199


def test_time_xor():
    one = vso.Time.dt((2010, 1, 1), (2010, 1, 2))
    a = one ^ vso.Time.dt((2010, 1, 1, 1), (2010, 1, 1, 2))
    
    assert a == attr.AttrOr(
        [vso.Time.dt((2010, 1, 1), (2010, 1, 1, 1)),
         vso.Time.dt((2010, 1, 1, 2), (2010, 1, 2))]
    )
    
    a ^= vso.Time.dt((2010, 1, 1, 4), (2010, 1, 1, 5))
    assert a == attr.AttrOr(
        [vso.Time.dt((2010, 1, 1), (2010, 1, 1, 1)),
         vso.Time.dt((2010, 1, 1, 2), (2010, 1, 1, 4)),
         vso.Time.dt((2010, 1, 1, 5), (2010, 1, 2))]
    )

def test_wave_xor():
    one = vso.Wave(0, 1000)
    a = one ^ vso.Wave(200, 400)
    
    assert a == attr.AttrOr([vso.Wave(0, 200), vso.Wave(400, 1000)])
    
    a ^= vso.Wave(600, 800)
    
    assert a == attr.AttrOr(
        [vso.Wave(0, 200), vso.Wave(400, 600), vso.Wave(800, 1000)])
