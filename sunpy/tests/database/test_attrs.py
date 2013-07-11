import pytest

from sunpy.database.database import Database
from sunpy.database import tables
from sunpy.database.attrs import Starred, Tag, walker
from sunpy.net.attr import DummyAttr, AttrAnd, AttrOr, and_


@pytest.fixture
def obj():
    return object()


@pytest.fixture
def session():
    database = Database('sqlite:///:memory:')
    database.create_tables()
    for i in xrange(1, 11):
        entry = tables.DatabaseEntry()
        database.add(entry)
        # every second entry gets starred
        if i % 2 == 0:
            database.star(entry)
        # every fifth entry gets the tag 'foo'
        if i % 5 == 0:
            database.tag(entry, 'foo')
    database.commit()
    return database.session


def test_starred_nonzero():
    assert Starred()


def test_starred_invert():
    assert not ~Starred()


def test_starred_and_different_types(obj):
    assert Starred() & obj == AttrAnd([Starred(), obj])


def test_starred_and_same_types():
    assert ~Starred() & ~Starred() == ~Starred()
    assert ~Starred() & Starred() == ~Starred()
    assert Starred() & ~Starred() == ~Starred()
    assert Starred() & Starred() == Starred()


def test_starred_or_different_types(obj):
    assert Starred() | obj == AttrOr([Starred(), obj])


def test_starred_or_same_types():
    assert ~Starred() | ~Starred() == ~Starred()
    assert ~Starred() | Starred() == Starred()
    assert Starred() | ~Starred() == Starred()
    assert Starred() | Starred() == Starred()


def test_starred_equality():
    assert Starred() == Starred()
    assert Starred() != ~Starred()
    assert ~Starred() != Starred()
    assert ~Starred() == ~Starred()


def test_starred_repr():
    assert repr(Starred()) == '<Starred(True)>'
    assert repr(~Starred()) == '<Starred(False)>'


def test_tag_repr():
    assert repr(Tag('foo')) == "<Tag('foo', False)>"


def test_walker_create_dummy(session):
    with pytest.raises(TypeError):
        walker.create(DummyAttr(), session)


def test_walker_create_starred_true(session):
    query = walker.create(Starred(), session)
    tag = tables.Tag('foo')
    tag.id = 1
    assert query.count() == 5
    assert query.all() == [
        tables.DatabaseEntry(id=2, starred=True),
        tables.DatabaseEntry(id=4, starred=True),
        tables.DatabaseEntry(id=6, starred=True),
        tables.DatabaseEntry(id=8, starred=True),
        tables.DatabaseEntry(id=10, starred=True, tags=[tag])]


def test_walker_create_starred_false(session):
    query = walker.create(~Starred(), session)
    tag = tables.Tag('foo')
    tag.id = 1
    assert query.count() == 5
    assert query.all() == [
        tables.DatabaseEntry(id=1),
        tables.DatabaseEntry(id=3),
        tables.DatabaseEntry(id=5, tags=[tag]),
        tables.DatabaseEntry(id=7),
        tables.DatabaseEntry(id=9)]


def test_walker_create_tag_positive(session):
    query = walker.create(Tag('foo'), session)
    tag = tables.Tag('foo')
    tag.id = 1
    assert query.count() == 2
    assert query.all() == [
        tables.DatabaseEntry(id=5, tags=[tag]),
        tables.DatabaseEntry(id=10, starred=True, tags=[tag])]


def test_walker_create_tag_negative(session):
    query = walker.create(~Tag('foo'), session)
    assert query.count() == 8
    assert query.all() == [
        tables.DatabaseEntry(id=1),
        tables.DatabaseEntry(id=2, starred=True),
        tables.DatabaseEntry(id=3),
        tables.DatabaseEntry(id=4, starred=True),
        tables.DatabaseEntry(id=6, starred=True),
        tables.DatabaseEntry(id=7),
        tables.DatabaseEntry(id=8, starred=True),
        tables.DatabaseEntry(id=9)]


@pytest.mark.xfail
def test_walker_create_anded_query(session):
    query = walker.create(and_(Tag('foo'), Starred()), session)
    assert query.count() == 1
    tag = tables.Tag('foo')
    tag.id = 1
    assert query.all() == [
        tables.DatabaseEntry(id=10, starred=True, tags=[tag])]
