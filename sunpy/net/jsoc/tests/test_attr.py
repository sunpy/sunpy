import pytest
import sunpy.net.jsoc as jsoc
import sunpy.net.jsoc.attrs as attrs
from sunpy.net.attr import Attr, AttrOr, AttrAnd


@pytest.mark.parametrize(("attr1, attr2"),
[  (attrs.Series('foo'), attrs.Series('boo')),
   (attrs.Protocol('a1'), attrs.Protocol('a2')),
   (attrs.Notify('email@somemail.com'), attrs.Notify('someemail@somemail.com')),
   (attrs.Compression('rice'), attrs.Compression('rice'))])
def test_and(attr1,attr2):
   pytest.raises(TypeError, lambda: attr1 & attr2)


def test_basicquery():
    a1 = attrs.Series('foo')
    t1 = attrs.Time('2012/01/01', '2013/1/2')
    ans1 = jsoc.jsoc.and_(a1, t1)
    assert isinstance(ans1, AttrAnd)
    assert len(ans1.attrs) == 2


def test_mediumquery():
    a1 = attrs.Series('foo1')
    a2 = attrs.Series('foo2')
    t1 = attrs.Time('2012/01/01', '2013/1/2')
    ans1 = jsoc.jsoc.and_(a1 | a2, t1)
    assert isinstance(ans1, AttrOr)
    assert isinstance(ans1.attrs[0], AttrAnd) 
    assert isinstance(ans1.attrs[1], AttrAnd)


def test_complexquery():
    a1 = attrs.Series('foo1')
    a2 = attrs.Series('foo2')
    t1 = attrs.Time('2012/01/01', '2013/1/2')
    t2 = attrs.Time('2012/01/01', '2013/1/3')
    ans1 = jsoc.jsoc.and_(a1 | a2, t1 | t2)
    assert isinstance(ans1.attrs[0], AttrOr)
    assert isinstance(ans1.attrs[0].attrs[0], AttrAnd)
    assert isinstance(ans1.attrs[0].attrs[1], AttrAnd)


