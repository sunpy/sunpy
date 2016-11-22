# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

from datetime import datetime

import pytest
from astropy import units as u

from sunpy.database.database import Database
from sunpy.database import tables
from sunpy.database.attrs import walker, Starred, Tag, Path, DownloadTime,\
    FitsHeaderEntry
from sunpy.net.attr import DummyAttr, AttrAnd, AttrOr
from sunpy.net import vso
from sunpy.extern.six.moves import range


@pytest.fixture
def obj():
    return object()


@pytest.fixture
def session():
    database = Database('sqlite:///:memory:')
    for i in range(1, 11):
        entry = tables.DatabaseEntry()
        database.add(entry)
        # every entry has a fake download time of 2005-06-15 i:00:00
        database.edit(entry, download_time=datetime(2005, 6, 15, i))
        # every second entry gets starred
        if i % 2 == 0:
            database.star(entry)
        # every third entry is stored in the path /tmp
        if i % 3 == 0:
            database.edit(entry, path='/tmp')
        # every fifth entry gets the tag 'foo'
        if i % 5 == 0:
            database.tag(entry, 'foo')
    # the last entry gets the FITS header entry INSTRUME=EIT
    entry.fits_header_entries.append(tables.FitsHeaderEntry('INSTRUME', 'EIT'))
    database.commit()
    return database.session


@pytest.fixture
def vso_session():
    client = vso.VSOClient()
    qr = client.query(
        vso.attrs.Time((2011, 9, 20, 1), (2011, 9, 20, 2)),
        vso.attrs.Instrument('RHESSI'))
    entries = tables.entries_from_query_result(qr)
    database = Database('sqlite:///:memory:')
    for entry in entries:
        database.add(entry)
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
    assert repr(Starred()) == '<Starred()>'
    assert repr(~Starred()) == '<~Starred()>'


def test_tag_repr():
    assert repr(Tag('foo')) == "<Tag('foo')>"
    assert repr(~Tag('foo')) == "<~Tag('foo')>"


def test_path_repr():
    assert repr(Path('/tmp')) == "<Path('/tmp')>"
    assert repr(~Path('/tmp')) == "<~Path('/tmp')>"


def test_downloadtime_repr():
    download_time = DownloadTime('2008-12-8', datetime(2009, 6, 12))
    expected_repr = (
        '<DownloadTime(datetime.datetime(2008, 12, 8, 0, 0), '
        'datetime.datetime(2009, 6, 12, 0, 0))>')
    assert repr(download_time) == expected_repr


def test_inverted_downloadtime_repr():
    download_time = ~DownloadTime('2008-12-8', datetime(2009, 6, 12))
    expected_repr = (
        '<~DownloadTime(datetime.datetime(2008, 12, 8, 0, 0), '
        'datetime.datetime(2009, 6, 12, 0, 0))>')
    assert repr(download_time) == expected_repr


def test_fitsheaderentry_repr():
    header_entry = FitsHeaderEntry('key', 'value')
    assert repr(header_entry) == "<FitsHeaderEntry('key', 'value')>"
    assert repr(~header_entry) == "<~FitsHeaderEntry('key', 'value')>"


def test_walker_create_dummy(session):
    with pytest.raises(TypeError):
        walker.create(DummyAttr(), session)


def test_walker_create_starred_true(session):
    entries = walker.create(Starred(), session)
    tag = tables.Tag('foo')
    tag.id = 1
    fits_header_entry = tables.FitsHeaderEntry('INSTRUME', 'EIT')
    fits_header_entry.id = 1
    assert len(entries) == 5
    assert entries == [
        tables.DatabaseEntry(
            id=2, starred=True, download_time=datetime(2005, 6, 15, 2)),
        tables.DatabaseEntry(
            id=4, starred=True, download_time=datetime(2005, 6, 15, 4)),
        tables.DatabaseEntry(
            id=6, path='/tmp', starred=True,
            download_time=datetime(2005, 6, 15, 6)),
        tables.DatabaseEntry(
            id=8, starred=True, download_time=datetime(2005, 6, 15, 8)),
        tables.DatabaseEntry(
            id=10, starred=True, tags=[tag],
            download_time=datetime(2005, 6, 15, 10),
            fits_header_entries=[fits_header_entry])]


def test_walker_create_starred_false(session):
    entries = walker.create(~Starred(), session)
    tag = tables.Tag('foo')
    tag.id = 1
    assert len(entries) == 5
    assert entries == [
        tables.DatabaseEntry(id=1, download_time=datetime(2005, 6, 15, 1)),
        tables.DatabaseEntry(
            id=3, path='/tmp', download_time=datetime(2005, 6, 15, 3)),
        tables.DatabaseEntry(
            id=5, tags=[tag], download_time=datetime(2005, 6, 15, 5)),
        tables.DatabaseEntry(
            id=7, download_time=datetime(2005, 6, 15, 7)),
        tables.DatabaseEntry(
            id=9, path='/tmp', download_time=datetime(2005, 6, 15, 9))]


def test_walker_create_tag_positive(session):
    entries = walker.create(Tag('foo'), session)
    tag = tables.Tag('foo')
    tag.id = 1
    fits_header_entry = tables.FitsHeaderEntry('INSTRUME', 'EIT')
    fits_header_entry.id = 1
    assert len(entries) == 2
    assert entries == [
        tables.DatabaseEntry(
            id=5, tags=[tag], download_time=datetime(2005, 6, 15, 5)),
        tables.DatabaseEntry(
            id=10, starred=True, tags=[tag],
            download_time=datetime(2005, 6, 15, 10),
            fits_header_entries=[fits_header_entry])]


def test_walker_create_tag_negative(session):
    entries = walker.create(~Tag('foo'), session)
    assert len(entries) == 8
    assert entries == [
        tables.DatabaseEntry(id=1, download_time=datetime(2005, 6, 15, 1)),
        tables.DatabaseEntry(
            id=2, starred=True, download_time=datetime(2005, 6, 15, 2)),
        tables.DatabaseEntry(
            id=3, path='/tmp', download_time=datetime(2005, 6, 15, 3)),
        tables.DatabaseEntry(
            id=4, starred=True, download_time=datetime(2005, 6, 15, 4)),
        tables.DatabaseEntry(
            id=6, path='/tmp', starred=True,
            download_time=datetime(2005, 6, 15, 6)),
        tables.DatabaseEntry(
            id=7, download_time=datetime(2005, 6, 15, 7)),
        tables.DatabaseEntry(
            id=8, starred=True, download_time=datetime(2005, 6, 15, 8)),
        tables.DatabaseEntry(
            id=9, path='/tmp', download_time=datetime(2005, 6, 15, 9))]


def test_walker_create_anded_query(session):
    entries = walker.create(Tag('foo') & Starred(), session)
    assert len(entries) == 1
    tag = tables.Tag('foo')
    tag.id = 1
    fits_header_entry = tables.FitsHeaderEntry('INSTRUME', 'EIT')
    fits_header_entry.id = 1
    assert tables.DatabaseEntry(
        id=10, starred=True, tags=[tag],
        download_time=datetime(2005, 6, 15, 10),
        fits_header_entries=[fits_header_entry]) in entries


def test_walker_create_ored_query(session):
    entries = walker.create(Tag('foo') | Starred(), session)
    assert len(entries) == 6
    tag = tables.Tag('foo')
    tag.id = 1
    fits_header_entry = tables.FitsHeaderEntry('INSTRUME', 'EIT')
    fits_header_entry.id = 1
    assert tables.DatabaseEntry(
        id=2, starred=True, download_time=datetime(2005, 6, 15, 2)) in entries
    assert tables.DatabaseEntry(
        id=4, starred=True, download_time=datetime(2005, 6, 15, 4)) in entries
    assert tables.DatabaseEntry(
        id=5, tags=[tag], download_time=datetime(2005, 6, 15, 5)) in entries
    assert tables.DatabaseEntry(
        id=6, path='/tmp', starred=True,
        download_time=datetime(2005, 6, 15, 6)) in entries
    assert tables.DatabaseEntry(
        id=8, starred=True, download_time=datetime(2005, 6, 15, 8)) in entries
    assert tables.DatabaseEntry(
        id=10, starred=True, tags=[tag],
        download_time=datetime(2005, 6, 15, 10),
        fits_header_entries=[fits_header_entry]) in entries


def test_walker_create_complex_query(session):
    query = Tag('foo') & Starred() | ~Tag('foo') & ~Starred()
    entries = walker.create(query, session)
    assert len(entries) == 5
    tag = tables.Tag('foo')
    tag.id = 1
    fits_header_entry = tables.FitsHeaderEntry('INSTRUME', 'EIT')
    fits_header_entry.id = 1
    assert tables.DatabaseEntry(
        id=1, download_time=datetime(2005, 6, 15, 1)) in entries
    assert tables.DatabaseEntry(
        id=3, path='/tmp', download_time=datetime(2005, 6, 15, 3)) in entries
    assert tables.DatabaseEntry(
        id=7, download_time=datetime(2005, 6, 15, 7)) in entries
    assert tables.DatabaseEntry(
        id=9, path='/tmp', download_time=datetime(2005, 6, 15, 9)) in entries
    assert tables.DatabaseEntry(
        id=10, starred=True, tags=[tag],
        download_time=datetime(2005, 6, 15, 10),
        fits_header_entries=[fits_header_entry]) in entries


def test_walker_create_path_attr_notfound(session):
    assert walker.create(Path('doesnotexist'), session) == []


def test_walker_create_path_attr_exists(session):
    entries = walker.create(Path('/tmp'), session)
    assert len(entries) == 3
    assert tables.DatabaseEntry(
        id=3, path='/tmp', download_time=datetime(2005, 6, 15, 3)) in entries
    assert tables.DatabaseEntry(
        id=6, path='/tmp', starred=True,
        download_time=datetime(2005, 6, 15, 6)) in entries
    assert tables.DatabaseEntry(
        id=9, path='/tmp', download_time=datetime(2005, 6, 15, 9)) in entries


def test_walker_create_path_inverted(session):
    tag = tables.Tag('foo')
    tag.id = 1
    fits_header_entry = tables.FitsHeaderEntry('INSTRUME', 'EIT')
    fits_header_entry.id = 1
    entries = walker.create(~Path('/tmp'), session)
    assert len(entries) == 7
    assert entries == [
        tables.DatabaseEntry(
            id=1, download_time=datetime(2005, 6, 15, 1)),
        tables.DatabaseEntry(
            id=2, starred=True, download_time=datetime(2005, 6, 15, 2)),
        tables.DatabaseEntry(
            id=4, starred=True, download_time=datetime(2005, 6, 15, 4)),
        tables.DatabaseEntry(
            id=5, tags=[tag], download_time=datetime(2005, 6, 15, 5)),
        tables.DatabaseEntry(id=7, download_time=datetime(2005, 6, 15, 7)),
        tables.DatabaseEntry(
            id=8, starred=True, download_time=datetime(2005, 6, 15, 8)),
        tables.DatabaseEntry(
            id=10, starred=True, tags=[tag],
            download_time=datetime(2005, 6, 15, 10),
            fits_header_entries=[fits_header_entry])]


def test_walker_create_downloadtime_notfound(session):
    download_time = DownloadTime(
        datetime(2005, 6, 15, 11), datetime(2005, 6, 15, 11))
    entries = walker.create(download_time, session)
    assert entries == []


def test_walker_create_downloadtime_exists(session):
    download_time = DownloadTime(
        datetime(2005, 6, 15, 7), datetime(2005, 6, 15, 9))
    entries = walker.create(download_time, session)
    assert entries == [
        tables.DatabaseEntry(id=7, download_time=datetime(2005, 6, 15, 7)),
        tables.DatabaseEntry(
            id=8, starred=True, download_time=datetime(2005, 6, 15, 8)),
        tables.DatabaseEntry(
            id=9, path='/tmp', download_time=datetime(2005, 6, 15, 9))]


def test_walker_create_downloadtime_inverted(session):
    tag = tables.Tag('foo')
    tag.id = 1
    fits_header_entry = tables.FitsHeaderEntry('INSTRUME', 'EIT')
    fits_header_entry.id = 1
    download_time = ~DownloadTime(
        datetime(2005, 6, 15, 7), datetime(2005, 6, 15, 9))
    entries = walker.create(download_time, session)
    assert len(entries) == 7
    assert entries == [
        tables.DatabaseEntry(
            id=1, download_time=datetime(2005, 6, 15, 1)),
        tables.DatabaseEntry(
            id=2, starred=True, download_time=datetime(2005, 6, 15, 2)),
        tables.DatabaseEntry(
            id=3, path='/tmp', download_time=datetime(2005, 6, 15, 3)),
        tables.DatabaseEntry(
            id=4, starred=True, download_time=datetime(2005, 6, 15, 4)),
        tables.DatabaseEntry(
            id=5, tags=[tag], download_time=datetime(2005, 6, 15, 5)),
        tables.DatabaseEntry(
            id=6, starred=True, path='/tmp',
            download_time=datetime(2005, 6, 15, 6)),
        tables.DatabaseEntry(
            id=10, starred=True, tags=[tag],
            download_time=datetime(2005, 6, 15, 10),
            fits_header_entries=[fits_header_entry])]


def test_walker_create_fitsheader(session):
    tag = tables.Tag('foo')
    tag.id = 1
    entries = walker.create(FitsHeaderEntry('INSTRUME', 'EIT'), session)
    fits_header_entry = tables.FitsHeaderEntry('INSTRUME', 'EIT')
    fits_header_entry.id = 1
    assert len(entries) == 1
    assert entries == [tables.DatabaseEntry(
        id=10, starred=True, tags=[tag],
        download_time=datetime(2005, 6, 15, 10),
        fits_header_entries=[fits_header_entry])]


def test_walker_create_fitsheader_inverted(session):
    tag = tables.Tag('foo')
    tag.id = 1
    entries = walker.create(~FitsHeaderEntry('INSTRUME', 'EIT'), session)
    assert len(entries) == 9
    assert entries == [
        tables.DatabaseEntry(
            id=1, download_time=datetime(2005, 6, 15, 1)),
        tables.DatabaseEntry(
            id=2, starred=True, download_time=datetime(2005, 6, 15, 2)),
        tables.DatabaseEntry(
            id=3, path='/tmp', download_time=datetime(2005, 6, 15, 3)),
        tables.DatabaseEntry(
            id=4, starred=True, download_time=datetime(2005, 6, 15, 4)),
        tables.DatabaseEntry(
            id=5, tags=[tag], download_time=datetime(2005, 6, 15, 5)),
        tables.DatabaseEntry(
            id=6, starred=True, path='/tmp',
            download_time=datetime(2005, 6, 15, 6)),
        tables.DatabaseEntry(id=7, download_time=datetime(2005, 6, 15, 7)),
        tables.DatabaseEntry(
            id=8, starred=True, download_time=datetime(2005, 6, 15, 8)),
        tables.DatabaseEntry(
            id=9, path='/tmp', download_time=datetime(2005, 6, 15, 9))]


@pytest.mark.online
def test_walker_create_vso_instrument(vso_session):
    entries = walker.create(vso.attrs.Instrument('RHESSI'), vso_session)
    assert entries == [
        tables.DatabaseEntry(id=1, source=u'RHESSI', provider=u'LSSP',
            physobs=u'intensity',
            fileid=u'/hessidata/2011/09/20/hsi_20110920_010920',
            observation_time_start=datetime(2011, 9, 20, 1, 9, 20),
            observation_time_end=datetime(2011, 9, 20, 2, 27, 40),
            instrument=u'RHESSI', size=-1.0, wavemin=0.4132806430668068,
            wavemax=7.293187818826002e-05),
        tables.DatabaseEntry(id=2, source=u'RHESSI', provider=u'LSSP',
            physobs=u'intensity',
            fileid=u'/hessidata/2011/09/19/hsi_20110919_233340',
            observation_time_start=datetime(2011, 9, 19, 23, 33, 40),
            observation_time_end=datetime(2011, 9, 20, 1, 9, 20),
            instrument=u'RHESSI', size=-1.0, wavemin=0.4132806430668068,
            wavemax=7.293187818826002e-05)]


@pytest.mark.online
def test_walker_create_wave(vso_session):
    entries = walker.create(vso.attrs.Wave(0 * u.AA, 10 * u.AA), vso_session)
    assert len(entries) == 2
    entries = walker.create(vso.attrs.Wave(5 * u.AA, 10 * u.AA), vso_session)
    assert len(entries) == 0


@pytest.mark.online
def test_walker_create_time(vso_session):
    time = vso.attrs.Time(
        datetime(2011, 9, 17, 0, 0, 0), datetime(2011, 9, 20, 0, 0, 0))
    entries = walker.create(time, vso_session)
    assert len(entries) == 1
    assert entries == [
        tables.DatabaseEntry(id=2, source=u'RHESSI', provider=u'LSSP',
            physobs=u'intensity',
            fileid=u'/hessidata/2011/09/19/hsi_20110919_233340',
            observation_time_start=datetime(2011, 9, 19, 23, 33, 40),
            observation_time_end=datetime(2011, 9, 20, 1, 9, 20),
            instrument=u'RHESSI', size=-1.0, wavemin=0.4132806430668068,
            wavemax=7.293187818826002e-05)]
