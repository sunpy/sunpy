import pytest

from sunpy.net import attr


class SA1(attr.SimpleAttr):
    pass


class SA2(attr.SimpleAttr):
    pass


class SA3(SA2):
    pass


class RA1(attr.Range):
    pass


@pytest.fixture
def walker():
    return attr.AttrWalker()


def test_creator(walker):
    CALLED = False

    @walker.add_creator(SA2)
    def new_creator(wlk, tree):
        nonlocal CALLED

        assert wlk is walker
        CALLED = True

    walker.create(SA2("Hello"))

    assert CALLED

    # Test that dispatch also works with subclasses
    CALLED = False

    walker.create(SA3("Hello"))

    assert CALLED


def test_applier(walker):
    CALLED = False

    @walker.add_applier(SA2)
    def new_applier(wlk, tree):
        nonlocal CALLED

        assert wlk is walker
        CALLED = True

    walker.apply(SA2("Hello"))

    assert CALLED

    CALLED = False

    # Test that dispatch also works with subclasses
    walker.apply(SA3("Hello"))

    assert CALLED


def test_creator_converter(walker):
    CALLED = False

    @walker.add_creator(SA1)
    def new_creator(wlk, tree):
        nonlocal CALLED

        assert wlk is walker
        CALLED = True

    with pytest.raises(TypeError):
        walker.create(RA1("Hello", "World"))

    assert not CALLED

    @walker.add_converter(RA1)
    def new_converter(arange):
        return SA1("CONVERTED")

    walker.create(RA1("Hello", "World"))

    assert CALLED


def test_applier_converter(walker):
    CALLED = False

    @walker.add_applier(SA1)
    def new_applier(wlk, tree):
        nonlocal CALLED

        assert wlk is walker
        CALLED = True

    with pytest.raises(TypeError):
        walker.apply(RA1("Hello", "World"))

    assert not CALLED

    @walker.add_converter(RA1)
    def new_converter(arange):
        return SA1("CONVERTED")

    walker.apply(RA1("Hello", "World"))

    assert CALLED
