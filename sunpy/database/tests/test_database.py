# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

import os
import glob
import shutil
import os.path
import configparser

import pytest
import sqlalchemy
from parfive.downloader import Downloader
from parfive.results import Results

from astropy import units
from astropy.utils.exceptions import AstropyUserWarning

import sunpy
import sunpy.data.test
from sunpy.data.test.waveunit import waveunitdir
from sunpy.database import (
    Database,
    EntryAlreadyAddedError,
    EntryAlreadyStarredError,
    EntryAlreadyUnstarredError,
    EntryNotFoundError,
    NoSuchTagError,
    TagAlreadyAssignedError,
    attrs,
    disable_undo,
    split_database,
    tables,
)
from sunpy.database.caching import LFUCache, LRUCache
from sunpy.database.commands import EmptyCommandStackError, NoSuchEntryError
from sunpy.database.tables import DatabaseEntry, FitsHeaderEntry, FitsKeyComment, JSONDump, Tag
from sunpy.io import fits
from sunpy.net import Fido
from sunpy.net import attrs as net_attrs
from sunpy.net import hek, vso

testpath = sunpy.data.test.rootdir
RHESSI_IMAGE = os.path.join(testpath, 'hsi_image_20101016_191218.fits')


"""
The 'hsi_image_20101016_191218.fits' file lies in the sunpy/data/test.
RHESSI_IMAGE  = sunpy/data/test/hsi_image_20101016_191218.fits

So, the tests in the database depends on the test under sunpy/data.
"""


@pytest.fixture
def database_using_lrucache():
    return Database('sqlite:///:memory:', LRUCache, cache_size=3)


@pytest.fixture
def database_using_lfucache():
    return Database('sqlite:///:memory:', LFUCache, cache_size=3)


@pytest.fixture
def database():
    return Database('sqlite:///:memory:')


@pytest.fixture
def fido_search_result():
    # A search query with responses from all instruments
    # No JSOC query
    return Fido.search(
        net_attrs.Time("2012/1/1", "2012/1/2"),
        net_attrs.Instrument('lyra') & net_attrs.Level.two | net_attrs.Instrument('eve') |
        net_attrs.Instrument('XRS') | net_attrs.Instrument('noaa-indices') |
        net_attrs.Instrument('noaa-predict') |
        (net_attrs.Instrument('norh') & net_attrs.Wavelength(17*units.GHz)) |
        (net_attrs.Instrument('rhessi') & net_attrs.Physobs("summary_lightcurve"))
    )


@pytest.fixture
def query_result():
    return vso.VSOClient().search(
        net_attrs.Time('20130801T200000', '20130801T200030'),
        net_attrs.Instrument('PLASTIC'))


@pytest.fixture
def download_qr():
    return vso.VSOClient().search(
        net_attrs.Time('2020-03-29', '2020-03-29'),
        net_attrs.Instrument('AIA'))


@pytest.fixture
def empty_query():
    return [
        net_attrs.Time((2012, 7, 3), (2012, 7, 4)),
        net_attrs.Instrument('EIT')]


@pytest.fixture
def download_query():
    return [
        net_attrs.Time((2013, 5, 19, 2), (2013, 5, 19, 2), (2013, 5, 19, 2)),
        net_attrs.Instrument('VIRGO') | net_attrs.Instrument('SECCHI')]


@pytest.fixture
def filled_database():
    database = Database('sqlite:///:memory:')
    for i in range(1, 11):
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


def test_config_url(monkeypatch):
    monkeypatch.setattr("sunpy.config", configparser.ConfigParser())
    url = 'sqlite:///'
    sunpy.config.add_section('database')
    sunpy.config.set('database', 'url', url)
    database = Database()
    assert database.url == url


def test_config_url_none(monkeypatch):
    monkeypatch.setattr("sunpy.config", configparser.ConfigParser())
    with pytest.raises(configparser.NoSectionError):
        Database()


def test_tags_unique(database):
    entry = DatabaseEntry()
    entry.tags = [Tag('foo')]
    database.add(entry)
    database.commit()
    entry.tags.append(Tag('foo'))
    with pytest.raises(sqlalchemy.orm.exc.FlushError):
        database.commit()


def test_setting_cache_size(database_using_lrucache):
    assert database_using_lrucache.cache_maxsize == 3
    assert database_using_lrucache.cache_size == 0
    for _ in range(5):
        database_using_lrucache.add(DatabaseEntry())
    assert len(database_using_lrucache) == 3
    assert database_using_lrucache.cache_size == 3
    assert database_using_lrucache.cache_maxsize == 3
    database_using_lrucache.set_cache_size(5)
    assert database_using_lrucache.cache_size == 3
    assert database_using_lrucache.cache_maxsize == 5
    for _ in range(5):
        database_using_lrucache.add(DatabaseEntry())
    assert len(database_using_lrucache) == 5
    assert database_using_lrucache.cache_size == 5
    assert database_using_lrucache.cache_maxsize == 5


def test_setting_cache_size_shrinking(database_using_lrucache):
    assert database_using_lrucache.cache_maxsize == 3
    assert database_using_lrucache.cache_size == 0
    for _ in range(5):
        database_using_lrucache.add(DatabaseEntry())
    assert len(database_using_lrucache) == 3
    assert database_using_lrucache.cache_maxsize == 3
    assert database_using_lrucache.cache_size == 3
    database_using_lrucache.set_cache_size(2)
    assert database_using_lrucache.cache_maxsize == 2
    assert database_using_lrucache.cache_size == 2
    assert len(database_using_lrucache) == 2
    assert list(database_using_lrucache) == [
        DatabaseEntry(id=4),
        DatabaseEntry(id=5)]
    for _ in range(5):
        database_using_lrucache.add(DatabaseEntry())
    assert len(database_using_lrucache) == 2
    assert database_using_lrucache.cache_maxsize == 2
    assert database_using_lrucache.cache_size == 2


def test_setting_cache_size_undo(database_using_lrucache):
    assert database_using_lrucache.cache_maxsize == 3
    assert database_using_lrucache.cache_size == 0
    for _ in range(5):
        database_using_lrucache.add(DatabaseEntry())
    assert len(database_using_lrucache) == 3
    database_using_lrucache.set_cache_size(1)
    assert database_using_lrucache.cache_size == 1
    assert len(database_using_lrucache) == 1
    database_using_lrucache.undo()
    assert len(database_using_lrucache) == 3


def test_get_entry_by_id_invalid(database):
    with pytest.raises(EntryNotFoundError):
        database.get_entry_by_id(0)


def test_get_entry_by_id_zero(filled_database):
    with pytest.raises(EntryNotFoundError):
        filled_database.get_entry_by_id(0)


def test_get_entry_by_id_accessible(filled_database):
    assert filled_database.get_entry_by_id(1) == DatabaseEntry(id=1)


def test_tags_property(database):
    assert database.tags == []


def test_get_existing_tag(database):
    entry = DatabaseEntry()
    database.tag(entry, 'tag')
    database.add(entry)
    expected_tag = Tag('tag')
    expected_tag.id = 1
    assert database.get_tag('tag') == expected_tag


def test_get_nonexting_tag(database):
    with pytest.raises(NoSuchTagError):
        database.get_tag('foo')


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


def test_tag_duplicate(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.tag(entry, 'tag')
    database.commit()
    with pytest.raises(TagAlreadyAssignedError):
        database.tag(entry, 'tag')


def test_tag_duplicates_before_adding(database):
    entry1 = DatabaseEntry()
    entry2 = DatabaseEntry()
    database.tag(entry1, 'tag')
    database.tag(entry2, 'tag')
    database.add(entry1)
    database.add(entry2)
    with pytest.raises(sqlalchemy.orm.exc.FlushError):
        database.commit()


def test_tag_undo(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.tag(entry, 'tag')
    assert len(database.tags) == 1
    assert len(entry.tags) == 1
    database.undo()
    assert len(entry.tags) == 0
    assert len(database.tags) == 0
    database.redo()
    assert len(database.tags) == 1
    assert len(entry.tags) == 1


def remove_nonexisting_tag(database):
    with pytest.raises(NoSuchTagError):
        database.remove_tag('foo')


def test_remove_tag(filled_database):
    foo = Tag('foo')
    foo.id = 1
    fourth_entry = filled_database.get_entry_by_id(4)
    assert foo in fourth_entry.tags
    eighth_entry = filled_database.get_entry_by_id(8)
    assert foo in eighth_entry.tags
    filled_database.remove_tag(fourth_entry, 'foo')
    assert foo not in fourth_entry.tags
    assert foo in filled_database.tags
    filled_database.remove_tag(eighth_entry, 'foo')
    assert foo not in eighth_entry.tags
    assert foo not in filled_database.tags


def test_remove_tag_undo_redo(filled_database):
    foo = Tag('foo')
    foo.id = 1
    fourth_entry = filled_database.get_entry_by_id(4)
    assert foo in fourth_entry.tags
    filled_database.remove_tag(fourth_entry, 'foo')
    assert foo not in fourth_entry.tags
    assert foo in filled_database.tags
    filled_database.undo()
    assert foo in fourth_entry.tags
    assert foo in filled_database.tags
    filled_database.redo()
    assert foo not in fourth_entry.tags
    assert foo in filled_database.tags
    eighth_entry = filled_database.get_entry_by_id(8)
    filled_database.remove_tag(eighth_entry, 'foo')
    assert foo not in eighth_entry.tags
    assert foo not in filled_database.tags
    filled_database.undo()
    assert foo not in fourth_entry.tags
    assert foo in eighth_entry.tags
    assert foo in filled_database.tags
    filled_database.redo()
    assert foo not in eighth_entry.tags
    assert foo not in filled_database.tags


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


def test_star_undo(database):
    entry = DatabaseEntry()
    assert not entry.starred
    database.star(entry)
    assert entry.starred
    database.undo()
    assert not entry.starred
    database.redo()
    assert entry.starred


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


def test_unstar_undo(database):
    entry = DatabaseEntry()
    entry.starred = True
    database.unstar(entry)
    assert not entry.starred
    database.undo()
    assert entry.starred
    database.redo()
    assert not entry.starred


def test_add_many(database):
    assert len(database) == 0
    database.add_many(DatabaseEntry() for _ in range(5))
    assert len(database) == 5
    database.undo()
    with pytest.raises(EmptyCommandStackError):
        database.undo()
    assert len(database) == 0
    database.redo()
    assert len(database) == 5


def test_add_many_with_existing_entry(database):
    evil_entry = DatabaseEntry()
    database.add(evil_entry)
    assert len(database) == 1
    with pytest.raises(EntryAlreadyAddedError):
        database.add_many([evil_entry])


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


@pytest.mark.remote_data
def test_add_entry_from_hek_qr(database):
    hek_res = hek.HEKClient().search(
        net_attrs.Time('2011/08/09 07:23:56', '2011/08/09 07:24:00'),
        hek.attrs.EventType('FL'))
    assert len(database) == 0
    database.add_from_hek_query_result(hek_res)
    # This number loves to change, so we are just going to test that it's added
    # *something*
    assert len(database) > 1


@pytest.mark.remote_data
def test_hek_query_download(monkeypatch, database, tmpdir):

    assert len(database) == 0

    records = ['94_1331820530-1331820530', '94_1331820542-1331820542',
               '94_1331820554-1331820554', '94_1331820566-1331820566',
               '94_1331820578-1331820578', '94_1331820590-1331820590',
               '94_1331820602-1331820602', '94_1331820614-1331820614',
               '94_1331820626-1331820626', '94_1331820638-1331820638']

    def mock_parfive_download(obj, *args, **kwargs):

        assert obj.queued_downloads == 10

        queue = obj.http_queue
        if not isinstance(queue, list):
            queue = list(queue._queue)
        obj_records = []

        for item in queue:
            url = item.keywords['url']
            obj_records.append(url[-24:])

        assert obj_records == records

        result = Results()
        result.append(str(tmpdir))
        return result

    def mock_entries_from_dir(*args, **kwargs):
        for i in range(10):
            yield DatabaseEntry()

    monkeypatch.setattr(Downloader, "download", mock_parfive_download)
    monkeypatch.setattr(tables, "entries_from_dir", mock_entries_from_dir)

    query = hek.HEKClient().search(
        net_attrs.Time('2019/03/10 14:40:10', '2019/04/11 16:40:50'),
        hek.attrs.EventType('FL')
    )

    database.download_from_hek_query_result(query[4], path=str(tmpdir))

    assert len(database) == 10


def num_entries_from_vso_query(db, query, path=None, file_pattern='',
                               overwrite=False):
    db.download_from_vso_query_result(
        query, path=path, overwrite=overwrite)
    fits_pattern = file_pattern
    num_of_fits_headers = sum(
        len(fits.get_header(file)) for file in glob.glob(fits_pattern))
    return num_of_fits_headers


@pytest.mark.remote_data
def test_vso_query_block_caching(database, download_qr, tmpdir):

    assert len(database) == 0

    # Download for all query response blocks and save the length
    # of database in num_of_fits_headers
    num_of_fits_headers = num_entries_from_vso_query(database, download_qr,
                                                     path=str(tmpdir.join('{file}.fits')),
                                                     file_pattern=str(tmpdir.join('*.fits')))

    assert len(database) == num_of_fits_headers and len(database) > 0

    # Emptying the database
    database.clear()
    database.commit()

    # Only downloading for the first query response block
    num_of_fits_headers_1 = num_entries_from_vso_query(database, download_qr[:1],
                                                       path=str(tmpdir.join('{file}.type1')),
                                                       file_pattern=str(tmpdir.join('*.type1')))

    assert len(database) == num_of_fits_headers_1 and len(database) > 0

    # Downloading for all query response blocks
    num_of_fits_headers_2 = num_entries_from_vso_query(database, download_qr,
                                                       path=str(tmpdir.join('{file}.type2')),
                                                       file_pattern=str(tmpdir.join('*.type2')))

    # Final length of the database should be the same as num_of_fits_headers.
    # This is done to ensure that the first query response block's files weren't
    # redownloaded. If they were redownloaded then length will be greater than
    # num_of_fits_headers as new entries are added to the database in case of a
    # download.

    assert len(database) == num_of_fits_headers_1 + num_of_fits_headers_2
    assert len(database) > 0

    assert num_of_fits_headers_1 + num_of_fits_headers_2 == num_of_fits_headers


@pytest.mark.remote_data
def test_vso_query_block_caching_with_overwrite_true_flag(database,
                                                          download_qr, tmpdir):

    assert len(database) == 0

    # Download for all query response blocks and save the length
    # of database in num_of_fits_headers

    num_of_fits_headers = num_entries_from_vso_query(database, download_qr,
                                                     path=str(tmpdir.join('{file}.fits')),
                                                     file_pattern=str(tmpdir.join('*.fits')))

    assert len(database) == num_of_fits_headers and len(database) > 0

    # Only downloading for the first query response block with caching disabled

    num_of_fits_headers_1 = num_entries_from_vso_query(database, download_qr[:1],
                                                       path=str(tmpdir.join('{file}.type1')),
                                                       file_pattern=str(tmpdir.join('*.type1')),
                                                       overwrite=True)

    # The files for the first query response block should be downloaded again
    # Old entries should be deleted, so len(database) should not change
    assert len(database) == num_of_fits_headers
    assert len(database) > 0


@pytest.mark.remote_data
def test_download_from_qr(database, download_qr, tmpdir):
    assert len(database) == 0
    database.download_from_vso_query_result(
        download_qr, path=str(tmpdir.join('{file}.fits')))
    fits_pattern = str(tmpdir.join('*.fits'))
    num_of_fits_headers = sum(
        len(fits.get_header(file)) for file in glob.glob(fits_pattern))
    assert len(database) == num_of_fits_headers > 0
    for entry in database:
        assert os.path.dirname(entry.path) == str(tmpdir)
    database.undo()
    assert len(database) == 0
    database.redo()
    assert len(database) == num_of_fits_headers > 0


@pytest.mark.remote_data
def test_add_entry_from_qr(database, query_result):
    assert len(database) == 0
    database.add_from_vso_query_result(query_result)
    assert len(database) == 16
    database.undo()
    assert len(database) == 0
    database.redo()
    assert len(database) == 16


@pytest.mark.remote_data
def test_add_entries_from_qr_duplicates(database, query_result):
    assert len(database) == 0
    database.add_from_vso_query_result(query_result)
    assert len(database) == 16
    with pytest.raises(EntryAlreadyAddedError):
        database.add_from_vso_query_result(query_result)


@pytest.mark.remote_data
def test_add_entries_from_qr_ignore_duplicates(database, query_result):
    assert len(database) == 0
    database.add_from_vso_query_result(query_result)
    assert len(database) == 16
    database.add_from_vso_query_result(query_result, True)
    assert len(database) == 32


@pytest.mark.remote_data
def test_add_entry_fido_search_result(database, fido_search_result):
    assert len(database) == 0
    database.add_from_fido_search_result(fido_search_result)
    assert len(database) == 66
    database.undo()
    assert len(database) == 0
    database.redo()
    assert len(database) == 66


@pytest.mark.remote_data
def test_add_entries_from_fido_search_result_JSOC_client(database):
    assert len(database) == 0
    search_result = Fido.search(
        net_attrs.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
        net_attrs.jsoc.Series('hmi.m_45s'),
        net_attrs.jsoc.Notify("sunpy@sunpy.org")
    )
    with pytest.raises(ValueError):
        database.add_from_fido_search_result(search_result)


@pytest.mark.remote_data
def test_add_entries_from_fido_search_result_duplicates(database, fido_search_result):
    assert len(database) == 0
    database.add_from_fido_search_result(fido_search_result)
    assert len(database) == 66
    with pytest.raises(EntryAlreadyAddedError):
        database.add_from_fido_search_result(fido_search_result)


@pytest.mark.remote_data
def test_add_entries_from_fido_search_result_ignore_duplicates(database, fido_search_result):
    assert len(database) == 0
    database.add_from_fido_search_result(fido_search_result)
    assert len(database) == 66
    database.add_from_fido_search_result(fido_search_result, True)
    assert len(database) == 2*66


def test_add_fom_path(database):
    assert len(database) == 0
    with pytest.warns(AstropyUserWarning, match='File may have been truncated'):
        database.add_from_dir(waveunitdir)
    assert len(database) == 4
    database.undo()
    assert len(database) == 0
    database.redo()
    assert len(database) == 4


def test_add_fom_path_duplicates(database):
    with pytest.warns(AstropyUserWarning, match='File may have been truncated'):
        database.add_from_dir(waveunitdir)
    assert len(database) == 4
    with pytest.raises(EntryAlreadyAddedError), pytest.warns(AstropyUserWarning, match='File may have been truncated'):
        database.add_from_dir(waveunitdir)


def test_add_fom_path_ignore_duplicates(database):
    with pytest.warns(AstropyUserWarning, match='File may have been truncated'):
        database.add_from_dir(waveunitdir)
    assert len(database) == 4
    with pytest.warns(AstropyUserWarning, match='File may have been truncated'):
        database.add_from_dir(waveunitdir, ignore_already_added=True)
    assert len(database) == 8


def test_add_from_file(database):
    assert len(database) == 0
    database.add_from_file(RHESSI_IMAGE)
    assert len(database) == 4
    # make sure that all entries have the same fileid
    fileid = database[0].fileid
    for entry in database:
        assert entry.fileid == fileid


def test_add_from_file_hdu_index(database):
    assert len(database) == 0
    database.add_from_file(RHESSI_IMAGE)
    assert len(database) == 4
    for i, entry in enumerate(database):
        assert entry.hdu_index == i


def test_add_from_file_duplicates(database):
    database.add_from_file(RHESSI_IMAGE)
    with pytest.raises(EntryAlreadyAddedError):
        database.add_from_file(RHESSI_IMAGE)


def test_add_from_file_ignore_duplicates(database):
    assert len(database) == 0
    database.add_from_file(RHESSI_IMAGE)
    assert len(database) == 4
    database.add_from_file(RHESSI_IMAGE, True)
    assert len(database) == 8


def test_edit_entry(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.commit()
    assert entry.id == 1
    database.edit(entry, id=42)
    assert entry.id == 42


def test_remove_many_entries(filled_database):
    bar = Tag('bar')
    bar.id = 2
    # required to check if `remove_many` adds any entries to undo-history
    filled_database.clear_histories()
    filled_database.remove_many(filled_database[:8])
    assert len(filled_database) == 2
    assert list(filled_database) == [
        DatabaseEntry(id=9),
        DatabaseEntry(id=10, tags=[bar])]
    filled_database.undo()
    assert len(filled_database) == 10
    with pytest.raises(EmptyCommandStackError):
        filled_database.undo()


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


def test_clear_empty_database(database):
    database.clear()


def test_clear_database(filled_database):
    assert len(filled_database) == 10
    filled_database.clear()
    assert not filled_database
    assert filled_database.session.query(JSONDump).all() == []
    assert filled_database.session.query(FitsHeaderEntry).all() == []
    assert filled_database.session.query(FitsKeyComment).all() == []
    assert filled_database.session.query(Tag).all() == []
    filled_database.undo()
    assert len(filled_database) == 10
    filled_database.redo()
    assert not filled_database
    assert filled_database.session.query(JSONDump).all() == []
    assert filled_database.session.query(FitsHeaderEntry).all() == []
    assert filled_database.session.query(FitsKeyComment).all() == []
    assert filled_database.session.query(Tag).all() == []


def test_getitem_notfound(database):
    with pytest.raises(IndexError):
        database[23]


def test_getitem_one(filled_database):
    first_entry = DatabaseEntry(id=1)
    assert filled_database[0] == first_entry


def test_getitem_getall(filled_database):
    entries = filled_database[:]
    assert entries == list(filled_database)


def test_getitem_custom(filled_database):
    entries = filled_database[1:5:2]
    foo = Tag('foo')
    foo.id = 1
    assert entries == [
        DatabaseEntry(id=2), DatabaseEntry(id=4, tags=[foo])]


def test_getitem_exceeding_range(filled_database):
    entries = filled_database[7:1000]
    foo = Tag('foo')
    foo.id = 1
    bar = Tag('bar')
    bar.id = 2
    assert entries == [
        DatabaseEntry(id=8, tags=[foo]),
        DatabaseEntry(id=9),
        DatabaseEntry(id=10, tags=[bar])]


def test_getitem_negative_index(filled_database):
    entry = filled_database[-4]
    assert entry == DatabaseEntry(id=7)


def test_getitem_negative_indices_slice(filled_database):
    entries = filled_database[-2:-8:-2]
    bar = Tag('bar')
    bar.id = 2
    assert entries == [
        DatabaseEntry(id=9),
        DatabaseEntry(id=7),
        DatabaseEntry(id=5, tags=[bar])]


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
    assert list(database_using_lrucache._cache.items()) == [
        (1, entry1), (2, entry2), (3, entry3)]
    database_using_lrucache.get_entry_by_id(1)
    database_using_lrucache.get_entry_by_id(3)
    entry4 = DatabaseEntry()
    database_using_lrucache.add(entry4)
    assert len(database_using_lrucache) == 3
    assert list(database_using_lrucache._cache.items()) == [
        (1, entry1), (3, entry3), (4, entry4)]


def test_lfu_cache(database_using_lfucache):
    assert not database_using_lfucache._cache
    entry1, entry2, entry3 = DatabaseEntry(), DatabaseEntry(), DatabaseEntry()
    database_using_lfucache.add(entry1)
    database_using_lfucache.add(entry2)
    database_using_lfucache.add(entry3)
    assert len(database_using_lfucache) == 3
    assert list(database_using_lfucache._cache.items()) == [
        (1, entry1), (2, entry2), (3, entry3)]
    # access the entries #1 and #2 to increment their counters
    database_using_lfucache.get_entry_by_id(1)
    database_using_lfucache.get_entry_by_id(2)
    entry4 = DatabaseEntry()
    database_using_lfucache.add(entry4)
    assert len(database_using_lfucache) == 3
    assert list(database_using_lfucache._cache.items()) == [
        (1, entry1), (2, entry2), (4, entry4)]


def test_query_missing_arg(database):
    with pytest.raises(TypeError):
        database.search()


def test_query_unexpected_kwarg(database):
    with pytest.raises(TypeError):
        database.search(attrs.Starred(), foo=42)


def test_query(filled_database):
    foo = Tag('foo')
    foo.id = 1
    bar = Tag('bar')
    bar.id = 2
    entries = filled_database.search(
        attrs.Tag('foo') | attrs.Tag('bar'), sortby='id')
    assert len(entries) == 4
    assert entries == [
        DatabaseEntry(id=4, tags=[foo]),
        DatabaseEntry(id=5, tags=[bar]),
        DatabaseEntry(id=8, tags=[foo]),
        DatabaseEntry(id=10, tags=[bar])]


def test_fetch_missing_arg(database):
    with pytest.raises(TypeError):
        database.fetch()


@pytest.mark.remote_data
def test_fetch_empty_query_result(database, empty_query):
    database.fetch(*empty_query)
    with pytest.raises(EmptyCommandStackError):
        database.undo()
    assert len(database) == 0


@pytest.mark.remote_data
def test_fetch(database, download_query, tmpdir):
    assert len(database) == 0
    database.default_waveunit = 'angstrom'
    database.fetch(
        *download_query, path=str(tmpdir.join('{file}.fits')), progress=True)
    fits_pattern = str(tmpdir.join('*.fits'))
    num_of_fits_headers = sum(
        len(fits.get_header(file)) for file in glob.glob(fits_pattern))
    assert len(database) == num_of_fits_headers
    for entry in database:
        assert os.path.dirname(entry.path) == str(tmpdir)
    database.undo()
    assert len(database) == 0
    database.redo()
    # Make this resilitent to vso changes while we chase this up with VSO 2018-03-07
    assert len(database) in (2, 4)


@pytest.mark.remote_data
def test_fetch_duplicates(database, download_query, tmpdir):
    assert len(database) == 0
    database.default_waveunit = 'angstrom'
    database.fetch(
        *download_query, path=str(tmpdir.join('{file}.fits')), progress=True)
    # FIXME The len(database) changes with time b/w 2 and 4.
    # Temp fix is  len(db) in (2, 4) until we find a better solution

    # 42 is the answer to life, the universe and everything ...
    # ... and this is no coincidence. So ...
    assert str(len(database)) in '42'

    download_time = database[0].download_time
    database.fetch(*download_query, path=str(tmpdir.join('{file}.fits')))
    # Make this resilitent to vso changes while we chase this up with VSO 2018-03-07
    assert len(database) in (2, 4)
    # The old file should be untouched because of the query result block
    # level caching
    assert database[0].download_time == download_time


def test_fetch_missing_arg(database):
    with pytest.raises(TypeError):
        database.fetch()


@pytest.mark.remote_data
def test_fetch(database, download_query, tmpdir):
    assert len(database) == 0
    database.default_waveunit = 'angstrom'
    database.fetch(*download_query, path=str(tmpdir.join('{file}.fits')))
    assert len(database) in (2, 4)
    download_time = database[0].download_time
    database.fetch(*download_query, path=str(tmpdir.join('{file}.fits')))
    assert len(database) in (2, 4)
    assert database[0].download_time == download_time


@pytest.mark.remote_data
def test_fetch_separate_filenames():
    # Setup
    db = Database('sqlite:///')

    download_query = [
        net_attrs.Time('2012-08-05', '2012-08-05 00:00:05'),
        net_attrs.Instrument('AIA')
    ]

    tmp_test_dir = os.path.join(
        sunpy.config.get('downloads', 'download_dir'),
        'tmp_test_dir/'
    )

    if not os.path.isdir(tmp_test_dir):
        os.makedirs(tmp_test_dir)

    path = tmp_test_dir + '{file}'

    db.fetch(*download_query, path=path)

    # Test
    assert len(db) == 4

    dir_contents = os.listdir(tmp_test_dir)
    assert 'aia_lev1_335a_2012_08_05t00_00_02_62z_image_lev1.fits' in dir_contents
    assert 'aia_lev1_94a_2012_08_05t00_00_01_12z_image_lev1.fits' in dir_contents
    assert os.path.isfile(os.path.join(
        tmp_test_dir, 'aia_lev1_335a_2012_08_05t00_00_02_62z_image_lev1.fits'))
    assert os.path.isfile(os.path.join(
        tmp_test_dir, 'aia_lev1_94a_2012_08_05t00_00_01_12z_image_lev1.fits'))

    # Teardown
    shutil.rmtree(tmp_test_dir)


@pytest.mark.remote_data
def test_disable_undo(database, download_query, tmpdir):
    entry = DatabaseEntry()
    with disable_undo(database) as db:
        db.set_cache_size(5)
        db.add(entry)
        db.commit()
        db.remove(entry)
        db.default_waveunit = 'angstrom'
        db.fetch(*download_query, path=str(tmpdir.join('{file}.fits')))
        entry = db[0]
        db.tag(entry, 'foo', 'bar')
        db.remove_tag(entry, 'foo')
        db.star(entry)
        db.unstar(entry)
        db.add_many([entry, entry], ignore_already_added=True)
        db.add(entry, ignore_already_added=True)
        db.add_from_dir(str(tmpdir))
        db.clear()
    with pytest.raises(EmptyCommandStackError):
        database.undo()


@pytest.fixture
def default_waveunit_database():
    unit_database = Database('sqlite:///:memory:',
                             default_waveunit=units.meter)
    str_database = Database('sqlite:///:memory:', default_waveunit="m")
    return unit_database, str_database


def test_default_waveunit(default_waveunit_database):
    unit_database, str_database = default_waveunit_database
    assert isinstance(unit_database.default_waveunit, units.UnitBase)
    assert isinstance(str_database.default_waveunit, units.UnitBase)


@pytest.fixture
def split_function_database():
    """
    Generates a custom database to test the split_database function
    """
    database = Database('sqlite:///:memory:')
    for i in range(1, 11):
        entry = DatabaseEntry()
        database.add(entry)
        # every fourth entry gets the instrument 'EIA'
        if i % 4 == 0:
            database.edit(entry, instrument='EIA')
        # every fifth entry gets the instrument 'AIA_3'
        elif i % 5 == 0:
            database.edit(entry, instrument='AIA_3')
        # every other entry gets instrument 'RHESSI'
        else:
            database.edit(entry, instrument='RHESSI')
        # all  entries have provider 'xyz'
        database.edit(entry, provider='xyz')
    database.commit()
    return database


def test_split_database(split_function_database, database):
    # Send all entries with instrument='EIA' to destination_database
    split_function_database, database = split_database(
        split_function_database, database, net_attrs.Instrument('EIA'))

    observed_source_entries = split_function_database.search(
        net_attrs.Provider('xyz'), sortby='id')
    observed_destination_entries = database.search(net_attrs.Provider('xyz'))

    assert observed_source_entries == [
        DatabaseEntry(id=1, instrument='RHESSI', provider='xyz'),
        DatabaseEntry(id=2, instrument='RHESSI', provider='xyz'),
        DatabaseEntry(id=3, instrument='RHESSI', provider='xyz'),
        DatabaseEntry(id=5, instrument='AIA_3', provider='xyz'),
        DatabaseEntry(id=6, instrument='RHESSI', provider='xyz'),
        DatabaseEntry(id=7, instrument='RHESSI', provider='xyz'),
        DatabaseEntry(id=9, instrument='RHESSI', provider='xyz'),
        DatabaseEntry(id=10, instrument='AIA_3', provider='xyz'),
    ]
    assert observed_destination_entries == [
        DatabaseEntry(id=4, instrument='EIA', provider='xyz'),
        DatabaseEntry(id=8, instrument='EIA', provider='xyz'),
    ]
