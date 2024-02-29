import pytest

import astropy.units as u

import sunpy.net.jsoc as jsoc
import sunpy.net.jsoc.attrs as attrs
from sunpy.net import _attrs as core_attrs
from sunpy.net.attr import AttrAnd, AttrOr


@pytest.mark.parametrize((("attr1", "attr2")),
                         [(attrs.Series('foo'), attrs.Series('boo')),
                          (attrs.Protocol('a1'), attrs.Protocol('a2')),
                          (attrs.Notify('email@somemail.com'),
                           attrs.Notify('someemail@somemail.com'))])
def test_and(attr1, attr2):
    with pytest.raises(TypeError, match=r"unsupported operand type\(s\) for \&"):
        attr1 & attr2


def test_basicquery():
    a1 = attrs.Series('foo')
    t1 = core_attrs.Time('2012/01/01', '2013/1/2')
    ans1 = jsoc.jsoc.and_(a1, t1)
    assert isinstance(ans1, AttrAnd)
    assert len(ans1.attrs) == 2


def test_mediumquery():
    a1 = attrs.Series('foo1')
    a2 = attrs.Series('foo2')
    t1 = core_attrs.Time('2012/01/01', '2013/1/2')
    ans1 = jsoc.jsoc.and_(a1 | a2, t1)
    assert isinstance(ans1, AttrOr)
    assert isinstance(ans1.attrs[0], AttrAnd)
    assert isinstance(ans1.attrs[1], AttrAnd)


def test_complexquery():
    a1 = attrs.Series('foo1')
    a2 = attrs.Series('foo2')
    t1 = core_attrs.Time('2012/01/01', '2013/1/2')
    t2 = core_attrs.Time('2012/01/01', '2013/1/3')
    ans1 = jsoc.jsoc.and_(a1 | a2, t1 | t2)
    assert isinstance(ans1.attrs[0], AttrOr)
    assert isinstance(ans1.attrs[0].attrs[0], AttrAnd)
    assert isinstance(ans1.attrs[0].attrs[1], AttrAnd)


def test_wavelength_error():
    with pytest.raises(TypeError):
        attrs.Wavelength('wobble')
    with pytest.raises(TypeError):
        attrs.Wavelength(3.24)
    with pytest.raises(TypeError):
        attrs.Wavelength((3, 3))


def test_wave_self():
    w1 = attrs.Wavelength(193*u.AA)
    assert jsoc.jsoc.and_(w1 | w1) is w1


def test_duplicate():
    w1 = attrs.Wavelength(193*u.AA)
    w2 = attrs.Wavelength(193*u.AA)
    assert jsoc.jsoc.and_(w1 | w2).min is w1.min


def test_random():
    w1 = attrs.Wavelength(193*u.AA)
    w2 = attrs.Series('spam')
    assert jsoc.jsoc.and_(w1 | w2) == AttrOr([w1, w2])


def test_empty_notify():
    with pytest.raises(ValueError, match="Notify attribute must contain an email address"):
        attrs.Notify(None)


def test_not_email_notify():
    with pytest.raises(ValueError, match="Notify attribute must contain an '@' symbol to be a valid email address"):
        attrs.Notify("someemailthatisntone")
