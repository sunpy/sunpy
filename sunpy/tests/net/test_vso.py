import pytest

from sunpy.net import vso

def pytest_funcarg__eit(request):
    return vso.Instrument('eit')


def pytest_funcarg__client(request):
    return vso.VSOClient()


def pytest_funcarg__iclient(request):
    return vso.InteractiveVSOClient()


def test_simpleattr_apply():
    attr = vso._SimpleAttr('test', 1)
    dct = {}
    attr.apply(dct)
    assert dct['test'] == 1


def test_simpleattr_create(client):
    attr = vso._SimpleAttr('instrument', 'eit')
    assert attr.create(client.api)[0].instrument == 'eit'


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
    tst = {'foo': 'a', 'bar': 'b'}
    attr = vso._ComplexAttr(['test'], tst)
    dct = {'test': {}}
    attr.apply(dct)
    assert dct['test'] == tst


def test_complexattr_create(client):
    attr = vso._ComplexAttr(['time'], {'start': 'test'})
    assert attr.create(client.api)[0].time.start == 'test'


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
    for name, factor in vso.Wave.energy:
        w = vso.Wave(62 / factor, 62 / factor, name)
        assert int(w.min) == 199
    
    w = vso.Wave(62, 62, 'eV')
    assert int(w.min) == 199
    w = vso.Wave(62e-3, 62e-3, 'keV')
    assert int(w.min) == 199

    for name, factor in vso.Wave.frequency:
        w = vso.Wave(1.506e16 / factor, 1.506e16 / factor, name)
        assert int(w.min) == 199
    
    w = vso.Wave(1.506e16, 1.506e16, 'Hz')
    assert int(w.min) == 199
    w = vso.Wave(1.506e7, 1.506e7, 'GHz')
    assert int(w.min) == 199
