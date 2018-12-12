# -*- coding: utf-8 -*-
from collections import defaultdict
import textwrap

import pytest

from sunpy.net import attr
from sunpy.net.attr import make_tuple

# TODO: Refactor the Attr tests, its too cluttered.

class Instrument(attr.SimpleAttr):
    """
    Dummy Instrument Class.
    """
    def __init__(self, value):
        super().__init__(value)


def EmptyAttr():
    attr.Attr._attr_registry = defaultdict(make_tuple)
    return attr.Attr


@pytest.fixture
def AIA():
    return Instrument('AIA')


@pytest.fixture
def NUM():
    return Instrument('1')


@pytest.fixture
def NUMBER():
    return Instrument('1AIA')


@pytest.fixture
def NUMBERS():
    return Instrument('12AIAs')


@pytest.fixture
def HMI():
    return Instrument('HMI')


@pytest.fixture
def SPEC():
    return Instrument('_!£!THIS_NAME!"!ISSPECIAL~~##')


@pytest.fixture
def KEYWORD():
    return Instrument('class')


class SA1(attr.SimpleAttr):
    pass


class SA2(attr.SimpleAttr):
    pass


class SA3(attr.SimpleAttr):
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


def test_attr_metamagic(AIA, HMI):
    attr.Attr.update_values({Instrument: [('AIA', 'This is AIA, it takes data')]})
    # .name is the attribute name return
    assert attr.Attr._attr_registry[Instrument].name == [AIA.value.lower()]
    # .name_long is the original name
    assert attr.Attr._attr_registry[Instrument].name_long == [AIA.value]
    # .des is the description of the item.
    assert attr.Attr._attr_registry[Instrument].desc == ['This is AIA, it takes data']
    # The _value_registry on the Attr object does not get cleaned.
    # So by adding it again to the same type, in this case Instrument the list is appended.
    attr.Attr.update_values({Instrument: [('HMI', 'This is HMI, it lives next to AIA')]})
    assert attr.Attr._attr_registry[Instrument].name == [AIA.value.lower(), HMI.value.lower()]
    assert attr.Attr._attr_registry[Instrument].name_long == [AIA.value, HMI.value]
    assert attr.Attr._attr_registry[Instrument].desc == ['This is AIA, it takes data', 'This is HMI, it lives next to AIA']

    # Tests the print out for the first two inputs only
    output = textwrap.dedent(
    "Attribute Name | Full Name | Description                      \n"
    "---------------+-----------+----------------------------------\n"
    "aia            | AIA       | This is AIA, it takes data       \n"
    "hmi            | HMI       | This is HMI, it lives next to AIA")

    assert str(Instrument) == output

    # Clean Registry
    EmptyAttr()


def test_attr_dynamic(AIA, HMI):
    # This checks the dynamic attribute creation.
    attr.Attr.update_values({Instrument: [('AIA', 'This is AIA, it takes data')]})
    attr.Attr.update_values({Instrument: [('HMI', 'This is HMI, it lives next to AIA')]})
    assert Instrument.aia == AIA
    assert Instrument.hmi == HMI
    # Clean Registry
    EmptyAttr()


def test_attr_dir(AIA, HMI):
    # Test for __dir__
    attr.Attr.update_values({Instrument: [('AIA', 'This is AIA, it takes data')]})
    attr.Attr.update_values({Instrument: [('HMI', 'This is HMI, it lives next to AIA')]})
    assert 'aia' in dir(Instrument)
    assert 'hmi' in dir(Instrument)
    # Clean Registry
    EmptyAttr()


def test_attr_sanity(SPEC):
    attr.Attr.update_values({Instrument: [('_!£!THIS_NAME!"!ISSPECIAL~~##', 'To test the attribute cleaning.')]})
    # This checks for sanitization of names.
    assert attr.Attr._attr_registry[Instrument].name == ['thisnameisspecial']
    assert attr.Attr._attr_registry[Instrument].name_long == ['_!£!THIS_NAME!"!ISSPECIAL~~##']
    assert attr.Attr._attr_registry[Instrument].desc == ['To test the attribute cleaning.']

    # Clean Registry
    EmptyAttr()


def test_attr_keyword(KEYWORD):
    attr.Attr.update_values({Instrument: [('class', 'Keyword checking.')]})
    # This checks for sanitization of names.
    assert attr.Attr._attr_registry[Instrument].name == ['class_']
    assert attr.Attr._attr_registry[Instrument].name_long == ['class']
    assert attr.Attr._attr_registry[Instrument].desc == ['Keyword checking.']

    # Clean Registry
    EmptyAttr()


def test_attr_num(NUM):
    attr.Attr.update_values({Instrument: [('1', 'One')]})
    # This checks for sanitization of names.
    assert attr.Attr._attr_registry[Instrument].name == ['one']
    assert attr.Attr._attr_registry[Instrument].name_long == ['1']
    assert attr.Attr._attr_registry[Instrument].desc == ['One']
    assert Instrument.one == NUM

    # Clean Registry
    EmptyAttr()


def test_attr_number(NUMBER):
    attr.Attr.update_values({Instrument: [('1AIA', 'One Number first.')]})
    # This checks for sanitization of names.
    assert attr.Attr._attr_registry[Instrument].name == ['one_aia']
    assert attr.Attr._attr_registry[Instrument].name_long == ['1AIA']
    assert attr.Attr._attr_registry[Instrument].desc == ['One Number first.']
    assert Instrument.one_aia == NUMBER

    # Clean Registry
    EmptyAttr()


def test_attr_numbes(NUMBERS):
    attr.Attr.update_values({Instrument: [('12AIAs', 'That is too many AIAs')]})
    # This checks for sanitization of names.
    assert attr.Attr._attr_registry[Instrument].name == ['12aias']
    assert attr.Attr._attr_registry[Instrument].name_long == ['12AIAs']
    assert attr.Attr._attr_registry[Instrument].desc == ['That is too many AIAs']
    assert '12aias' in dir(Instrument)

    # Clean Registry
    EmptyAttr()


def test_attr_iterable_length(AIA):
    # not iterable
    with pytest.raises(ValueError):
        attr.Attr.update_values({Instrument: 'AIA'})
    # too many items
    with pytest.raises(ValueError):
        attr.Attr.update_values({Instrument: [('AIA', 'AIA is Nice', 'Error now')]})

    # Clean Registry
    EmptyAttr()
