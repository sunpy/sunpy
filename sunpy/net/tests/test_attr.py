from collections import defaultdict

import pytest

from sunpy.net import attr
from sunpy.net.attr import AttrMeta, make_tuple
from sunpy.net.dataretriever import GenericClient


class Instrument(attr.SimpleAttr):
    """
    Dummy Instrument Class.
    """


class Time(attr.Range):
    """
    Dummy Time Class.
    """


def EmptyAttr():
    AttrMeta._attr_registry = defaultdict(make_tuple)


@pytest.fixture
def ALL():
    return Instrument('all')


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
def POINTNUMBER():
    return Instrument('1.5')


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


class SA4(attr.SimpleAttr):
    pass


def test_empty():
    class TestAttr(attr.Attr):
        pass

    assert repr(TestAttr)


@pytest.mark.parametrize("different_type", [
    int, str, float, list, set, tuple, dict, object
])
def test_empty(different_type):
    attr_ = attr.Attr()
    assert attr_ != different_type()


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


def test_attr_multi_and_AttrAnd():
    a1 = SA1(1)
    a2 = SA2(2)
    a3 = SA3(3)
    a4 = SA4(4)
    a_and1 = (a2 & a3)
    a_and2 = (a1 & a4)
    an = a_and1 & a_and2

    assert isinstance(a_and1, attr.AttrAnd)
    assert isinstance(a_and2, attr.AttrAnd)
    assert isinstance(an, attr.AttrAnd)
    assert a1 in an.attrs
    assert a2 in an.attrs
    assert a3 in an.attrs
    assert a4 in an.attrs
    assert len(an.attrs) == 4


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
    # {cls: cls.register_values()}
    attr.Attr.update_values({GenericClient: {Instrument: [('AIA', 'This is AIA, it takes data')]}})
    # .name is the attribute name return
    assert attr.Attr._attr_registry[Instrument].name == [AIA.value.lower()]
    # .name_long is the original name
    assert attr.Attr._attr_registry[Instrument].name_long == [AIA.value]
    # .des is the description of the item.
    assert attr.Attr._attr_registry[Instrument].desc == ['This is AIA, it takes data']
    # The _value_registry on the Attr object does not get cleaned.
    # So by adding it again to the same type, in this case Instrument the list is appended.
    attr.Attr.update_values(
        {GenericClient: {Instrument: [('HMI', 'This is HMI, it lives next to AIA')]}})
    assert attr.Attr._attr_registry[Instrument].name == [AIA.value.lower(), HMI.value.lower()]
    assert attr.Attr._attr_registry[Instrument].name_long == [AIA.value, HMI.value]
    assert attr.Attr._attr_registry[Instrument].desc == [
        'This is AIA, it takes data', 'This is HMI, it lives next to AIA']

    # Tests the print out for the first two inputs only
    output = 'sunpy.net.tests.test_attr.Instrument\n\nDummy Instrument Class.\n\n\nAttribute Name  Client Full Name            Description           \n-------------- ------- --------- ---------------------------------\naia            Generic AIA       This is AIA, it takes data       \nhmi            Generic HMI       This is HMI, it lives next to AIA'
    assert str(Instrument) == output


def test_attr_dynamic(AIA, HMI):
    # This checks the dynamic attribute creation.
    attr.Attr.update_values({GenericClient: {Instrument: [('AIA', 'This is AIA, it takes data')]}})
    attr.Attr.update_values(
        {GenericClient: {Instrument: [('HMI', 'This is HMI, it lives next to AIA')]}})
    assert Instrument.aia == AIA
    assert Instrument.hmi == HMI


def test_attr_dir():
    # Test for __dir__
    attr.Attr.update_values({GenericClient: {Instrument: [('AIA', 'This is AIA, it takes data')]}})
    attr.Attr.update_values(
        {GenericClient: {Instrument: [('HMI', 'This is HMI, it lives next to AIA')]}})

    assert 'aia' in dir(Instrument)
    assert 'hmi' in dir(Instrument)


def test_attr_sanity():
    attr.Attr.update_values(
        {GenericClient: {Instrument: [('_!£!THIS_NAME!"!ISSPECIAL~~##', 'To test the attribute cleaning.')]}})
    # This checks for sanitization of names.
    assert '___this_name___isspecial____' in attr.Attr._attr_registry[Instrument].name
    assert '_!£!THIS_NAME!"!ISSPECIAL~~##' in attr.Attr._attr_registry[Instrument].name_long
    assert 'To test the attribute cleaning.' in attr.Attr._attr_registry[Instrument].desc


def test_attr_keyword():
    attr.Attr.update_values({GenericClient: {Instrument: [('class', 'Keyword checking.')]}})
    # This checks for sanitization of names.
    assert 'class_' in attr.Attr._attr_registry[Instrument].name
    assert 'class' in attr.Attr._attr_registry[Instrument].name_long
    assert 'Keyword checking.' in attr.Attr._attr_registry[Instrument].desc


def test_attr_num(NUM):
    attr.Attr.update_values({GenericClient: {Instrument: [('1', 'One')]}})
    # This checks for sanitization of names.
    assert 'one' in attr.Attr._attr_registry[Instrument].name
    assert '1' in attr.Attr._attr_registry[Instrument].name_long
    assert 'One' in attr.Attr._attr_registry[Instrument].desc
    assert Instrument.one == NUM


def test_attr_number(NUMBER):
    attr.Attr.update_values({GenericClient: {Instrument: [('1AIA', 'One Number first.')]}})
    # This checks for sanitization of names.
    assert 'one_aia' in attr.Attr._attr_registry[Instrument].name
    assert '1AIA' in attr.Attr._attr_registry[Instrument].name_long
    assert 'One Number first.' in attr.Attr._attr_registry[Instrument].desc
    assert Instrument.one_aia == NUMBER


def test_attr_number_point(POINTNUMBER):
    attr.Attr.update_values({GenericClient: {Instrument: [('1.5', 'One Point Five.')]}})
    # This checks for sanitization of names.
    assert 'onepointfive' in attr.Attr._attr_registry[Instrument].name
    assert '1.5' in attr.Attr._attr_registry[Instrument].name_long
    assert 'One Point Five.' in attr.Attr._attr_registry[Instrument].desc
    assert Instrument.onepointfive == POINTNUMBER


def test_attr_numbes():
    attr.Attr.update_values({GenericClient: {Instrument: [('12AIAs', 'That is too many AIAs')]}})
    # This checks for sanitization of names.
    assert 'one_2aias' in attr.Attr._attr_registry[Instrument].name
    assert '12AIAs' in attr.Attr._attr_registry[Instrument].name_long
    assert 'That is too many AIAs' in attr.Attr._attr_registry[Instrument].desc
    assert 'one_2aias' in dir(Instrument)


def test_attr_iterable_length():
    # not iterable
    with pytest.raises(ValueError):
        attr.Attr.update_values({GenericClient: {Instrument: 'AIA'}})
    # too many items
    with pytest.raises(ValueError):
        attr.Attr.update_values(
            {GenericClient: {Instrument: [('AIA', 'AIA is Nice', 'Error now')]}})


def test_asterisk_attrs(ALL):
    # This checks we can submit * to mean all attrs.
    attr.Attr.update_values({GenericClient: {Instrument: [('*')]}})
    assert Instrument.all == ALL
    assert "Instrument(all: All values of this type are supported.)" in repr(Instrument.all)


@pytest.mark.parametrize("wrong_name", [
    ("not star",), ("*whoops",)
])
def test_single_pair_argument_attrs(wrong_name):
    # This checks that other single string entries fail.
    with pytest.raises(ValueError):
        attr.Attr.update_values({GenericClient: {Instrument: [wrong_name]}})


def test_asterisk_attrs_time():
    # This checks we can submit * for time/wavelength (both are ranges)
    attr.Attr.update_values({GenericClient: {Time: [('*')]}})
    assert "all       All values of this type are supported." in repr(Time)
