from __future__ import absolute_import

import pytest

from sunpy.database import Database, EntryAlreadyAddedError,\
    EntryAlreadyStarredError, EntryAlreadyUnstarredError
from sunpy.database.tables import DatabaseEntry, Tag
from sunpy.database.commands import NoSuchEntryError
from sunpy.database.caching import LRUCache, LFUCache
from sunpy.database import attrs


@pytest.fixture
def database_without_tables():
    return Database('sqlite:///:memory:')


@pytest.fixture
def database_using_lrucache():
    d = Database('sqlite:///:memory:', LRUCache, cache_size=3)
    d.create_tables()
    return d


@pytest.fixture
def database_using_lfucache():
    d = Database('sqlite:///:memory:', LFUCache, cache_size=3)
    d.create_tables()
    return d


@pytest.fixture
def database():
    d = Database('sqlite:///:memory:')
    d.create_tables()
    return d


@pytest.fixture
def filled_database():
    database = Database('sqlite:///:memory:')
    database.create_tables()
    for i in xrange(1, 11):
        entry = DatabaseEntry()
        database.add(entry)
        # every fourth entry gets the tag 'foo'
        if i % 4 == 0:
            database.tag(entry, 'foo')
        # every fifth entry gets the tag 'bar'
        if i % 5 == 0:
            database.tag(entry, 'bar')
    database.commit()
    return database


def test_create_tables(database_without_tables):
    assert not database_without_tables._engine.has_table('data')
    database_without_tables.create_tables()
    assert database_without_tables._engine.has_table('data')


def test_tags_property(database):
    assert database.tags == []


def test_get_existing_tag(database):
    entry = DatabaseEntry()
    database.tag(entry, 'tag')
    database.add(entry)
    # dirty trick:
    expected_tag = Tag('tag')
    expected_tag.id = 1
    assert database.get_tag('tag') == expected_tag


def test_get_nonexting_tag(database):
    assert database.get_tag('foo') is None


def test_tag_missing_tags_arg(database):
    with pytest.raises(TypeError):
        database.tag(DatabaseEntry())


def test_tag_new_tag(database):
    entry = DatabaseEntry()
    database.tag(entry, 'tag')
    assert len(entry.tags) == 1
    database.add(entry)
    assert len(database.tags) == 1
    tag = entry.tags[0]
    assert tag.name == 'tag'
    assert tag in database.tags


def test_tag_existing_tag(database):
    entry1 = DatabaseEntry()
    entry2 = DatabaseEntry()
    database.tag(entry1, 'tag')
    database.tag(entry2, 'tag')
    assert entry1.tags == entry2.tags


# the following test raises a sqlalchemy.exc.IntegrityError exception because
# the uniqueness of the tag names can only be checked if the according entries
# are already saved in the database (otherwise the tags have no IDs)
@pytest.mark.xfail
def test_tag_duplicates_before_adding(database):
    entry1 = DatabaseEntry()
    entry2 = DatabaseEntry()
    database.tag(entry1, 'tag')
    database.tag(entry2, 'tag')
    database.add(entry1)
    database.add(entry2)
    database.commit()


def test_get_by_tags_missing_tags_arg(database):
    with pytest.raises(TypeError):
        database.get_by_tags()


def test_get_by_tags_no_matches(database):
    assert database.get_by_tags('foobar') == []


def test_get_by_tags_matching(database):
    entry1 = DatabaseEntry()
    entry2 = DatabaseEntry()
    database.add(entry1)
    database.add(entry2)
    database.tag(entry1, 'one')
    database.tag(entry2, 'one', 'two')
    assert database.get_by_tags('two') == [entry2]
    assert database.get_by_tags('one') == [entry1, entry2]


def test_star_entry(database):
    entry = DatabaseEntry()
    assert not entry.starred
    database.star(entry)
    assert entry.starred


def test_star_already_starred_entry(database):
    entry = DatabaseEntry()
    database.star(entry)
    with pytest.raises(EntryAlreadyStarredError):
        database.star(entry)


def unstar_entry(database):
    entry = DatabaseEntry()
    assert not entry.starred
    database.star(entry)
    assert entry.starred
    database.unstar(entry)
    assert not entry.starred


def test_unstar_already_unstarred_entry(database):
    with pytest.raises(EntryAlreadyUnstarredError):
        database.unstar(DatabaseEntry())


def test_unstar_already_unstarred_entry_ignore(database):
    entry = DatabaseEntry()
    database.unstar(entry, True)
    assert not entry.starred


def test_get_starred_entries(database):
    entry1 = DatabaseEntry()
    entry2 = DatabaseEntry()
    database.add(entry1)
    database.add(entry2)
    assert list(database.get_starred()) == []
    database.star(entry2)
    assert list(database.get_starred()) == [entry2]


def test_add_entry(database):
    entry = DatabaseEntry()
    assert entry.id is None
    database.add(entry)
    database.commit()
    assert entry.id == 1


def test_add_already_existing_entry(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.commit()
    with pytest.raises(EntryAlreadyAddedError):
        database.add(entry)


def test_add_already_existing_entry_ignore(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.add(entry, True)
    database.commit()
    assert entry.id == 1


def test_edit_entry(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.commit()
    assert entry.id == 1
    database.edit(entry, id=42)
    assert entry.id == 42


def test_remove_existing_entry(database):
    entry = DatabaseEntry()
    database.add(entry)
    assert database.session.query(DatabaseEntry).count() == 1
    assert entry.id == 1
    database.remove(entry)
    assert database.session.query(DatabaseEntry).count() == 0


def test_remove_nonexisting_entry(database):
    with pytest.raises(NoSuchEntryError):
        database.remove(DatabaseEntry())


def test_contains_exists(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.commit()
    assert entry in database


def test_contains_precommit(database):
    entry = DatabaseEntry()
    database.add(entry)
    assert entry not in database


def test_contains_notexists(database):
    assert DatabaseEntry() not in database


def test_iter(database):
    entry1 = DatabaseEntry()
    entry2 = DatabaseEntry()
    database.add(entry1)
    database.add(entry2)
    expected_entries = [entry1, entry2]
    entries = list(database)
    assert entries == expected_entries


def test_len(database):
    assert len(database) == 0
    database.session.add(DatabaseEntry())
    assert len(database) == 1


def test_lru_cache(database_using_lrucache):
    assert not database_using_lrucache._cache
    entry1, entry2, entry3 = DatabaseEntry(), DatabaseEntry(), DatabaseEntry()
    database_using_lrucache.add(entry1)
    database_using_lrucache.add(entry2)
    database_using_lrucache.add(entry3)
    assert len(database_using_lrucache) == 3
    entries = list(database_using_lrucache)
    assert entries[0] == entry1
    assert entries[1] == entry2
    assert entries[2] == entry3
    #assert database_using_lrucache._cache.items() == [
    #    (1, entry1), (2, entry2), (3, entry3)]
    database_using_lrucache.get_entry_by_id(1)
    database_using_lrucache.get_entry_by_id(3)
    entry4 = DatabaseEntry()
    database_using_lrucache.add(entry4)
    assert len(database_using_lrucache) == 3
    entries = list(database_using_lrucache)
    assert entries[0] == entry1
    assert entries[1] == entry3
    assert entries[2] == entry4
    #assert database_using_lrucache._cache.items() == [
    #    (1, entry1), (3, entry3), (4, entry4)]


def test_lfu_cache(database_using_lfucache):
    assert not database_using_lfucache._cache
    entry1, entry2, entry3 = DatabaseEntry(), DatabaseEntry(), DatabaseEntry()
    database_using_lfucache.add(entry1)
    database_using_lfucache.add(entry2)
    database_using_lfucache.add(entry3)
    assert len(database_using_lfucache) == 3
    entries = list(database_using_lfucache)
    assert entries[0] == entry1
    assert entries[1] == entry2
    assert entries[2] == entry3
    #assert database_using_lrucache._cache.items() == [
    #    (1, entry1), (2, entry2), (3, entry3)]
    # access the entries #1 and #2 to increment their counters
    database_using_lfucache.get_entry_by_id(1)
    database_using_lfucache.get_entry_by_id(2)
    entry4 = DatabaseEntry()
    database_using_lfucache.add(entry4)
    assert len(database_using_lfucache) == 3
    entries = list(database_using_lfucache)
    assert entries[0] == entry1
    assert entries[1] == entry2
    assert entries[2] == entry4
    #assert database_using_lrucache._cache.items() == [
    #    (1, entry1), (2, entry2), (4, entry4)]


def test_query_missing_arg(database):
    with pytest.raises(TypeError):
        database.query()


def test_query_unexpected_kwarg(database):
    with pytest.raises(TypeError):
        database.query(attrs.Starred(), foo=42)


def test_query(filled_database):
    foo = Tag('foo')
    foo.id = 1
    bar = Tag('bar')
    bar.id = 2
    entries = filled_database.query(attrs.Tag('foo') | attrs.Tag('bar'), sortby='id')
    assert len(entries) == 4
    assert entries == [
        DatabaseEntry(id=4, tags=[foo]),
        DatabaseEntry(id=5, tags=[bar]),
        DatabaseEntry(id=8, tags=[foo]),
        DatabaseEntry(id=10, tags=[bar])]
