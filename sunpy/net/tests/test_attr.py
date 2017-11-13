# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pytest

from sunpy.net import attr


class SA1(attr.SimpleAttr):
    pass


class SA2(attr.SimpleAttr):
    pass


class SA3(attr.SimpleAttr):
    pass


class SA4(attr.SimpleAttr):
    pass


def test_attr_and():
    a1 = SA1(1)
    a2 = SA2(2)

    an = a1 & a2

    assert isinstance(an, attr.AttrAnd)
    assert a1 in an.attrs
    assert a2 in an.attrs
    assert len(an.attrs) == 2


def test_attr_and_AttrAnd():
    a1 = SA1(1)
    a2 = SA2(2)
    a3 = SA3(3)

    an = a1 & (a2 & a3)

    assert isinstance(an, attr.AttrAnd)
    assert a1 in an.attrs
    assert a2 in an.attrs
    assert a3 in an.attrs
    assert len(an.attrs) == 3


def test_attr_and_AttrOr():
    a1 = SA1(1)
    a2 = SA2(2)
    a3 = SA3(3)

    an = a1 & (a2 | a3)

    assert isinstance(an, attr.AttrOr)
    for a in an.attrs:
        assert isinstance(a, attr.AttrAnd)
    assert len(an.attrs) == 2


def test_attr_hash():
    a1 = SA1(1)
    a2 = SA1(1)
    a3 = SA1(3)

    assert hash(a1) == hash(a2)
    assert hash(a3) != hash(a1)


def test_attr_collies():
    a1 = attr.Attr()
    with pytest.raises(NotImplementedError):
        a1.collides(1)


def test_attr_or():
    a1 = SA1(1)
    a2 = SA2(2)

    an = a1 | a2

    assert isinstance(an, attr.AttrOr)
    assert a1 in an.attrs
    assert a2 in an.attrs
    assert len(an.attrs) == 2

    a1 = SA1(1)
    a2 = SA2(1)

    an = a1 | a2

    assert an is a1


def test_simpleattr_collides():
    a1 = SA1(1)

    with pytest.raises(TypeError):
        a1 & a1


def test_simple_attr_repr():
    a1 = SA1("test string")

    assert "test string" in repr(a1)
    assert "SA1" in repr(a1)


def test_dummyattr():
    one = attr.DummyAttr()
    other = attr.ValueAttr({'a': 'b'})
    assert (one | other) is other
    assert (one & other) is other


def test_dummyattr_hash():
    one = attr.DummyAttr()
    assert hash(one) == hash(None)


def test_dummyattr_collides():
    one = attr.DummyAttr()
    two = attr.DummyAttr()
    assert one.collides(two) is False


def test_dummyattr_eq():
    one = attr.DummyAttr()
    two = attr.DummyAttr()
    other = attr.ValueAttr({'a': 'b'})
    assert one == two
    assert one != other


def test_and_nesting():
    a1 = SA1(1)
    a2 = SA2(2)
    a3 = SA3(3)

    a = attr.and_(a1, attr.AttrAnd((a2, a3)))
    # Test that the nesting has been removed.
    assert len(a.attrs) == 3


def test_or_nesting():
    a1 = SA1(1)
    a2 = SA2(2)
    a3 = SA3(3)

    a = attr.or_(a1, attr.AttrOr((a2, a3)))
    # Test that the nesting has been removed.
    assert len(a.attrs) == 3


def test_attrand_repr():
    a1 = SA1(1)
    a2 = SA2(2)
    aand = a1 & a2

    assert "SA1" in repr(aand)
    assert "SA2" in repr(aand)


def test_attrand_eq():
    a1 = SA1(1)
    a2 = SA2(2)
    aand = a1 & a2
    aand2 = a1 & a2

    assert not aand == a1
    assert aand == aand2


def test_attrand_hash():
    a1 = SA1(1)
    a2 = SA2(2)
    a3 = SA2(3)
    aand = a1 & a2
    aand2 = a1 & a2
    aand3 = a1 & a3

    assert hash(aand) == hash(aand2)
    assert hash(aand) != hash(aand3)


def test_attrand_collides():
    a1 = SA1(1)
    a2 = SA2(2)
    a3 = SA2(3)
    aand = a1 & a2
    aand2 = a1 & a2

    assert aand.collides(a3)
    assert not aand.collides(aand2)


def test_attrand_and():
    a1 = SA1(1)
    a2 = SA2(2)
    a3 = SA3(3)
    a4 = SA4(4)
    aand = a1 & a2
    aand2 = a4 & a3
    aor = a4 | a3

    with pytest.raises(TypeError):
        aand & aand

    assert isinstance(aand & aand2, attr.AttrAnd)

    assert isinstance(aand & aor, attr.AttrOr)

    assert isinstance(aand & a3, attr.AttrAnd)
    assert len((aand & a3).attrs) == 3


def test_attror_or():
    a1 = SA1(1)
    a2 = SA2(2)
    a3 = SA3(3)
    a4 = SA4(4)
    aor = a1 | a2
    aor2 = a4 | a3

    comp = aor | aor2
    assert isinstance(comp, attr.AttrOr)
    assert len(comp.attrs) == 4

    assert len((aor | a3).attrs) == 3


def test_attror_and():
    a1 = SA1(1)
    a2 = SA2(2)
    a3 = SA3(3)
    aor = a1 | a2

    comp = aor & a3

    assert isinstance(comp, attr.AttrOr)
    assert len(comp.attrs) == 2
    assert isinstance(comp.attrs[0], attr.AttrAnd)
    assert isinstance(comp.attrs[1], attr.AttrAnd)
