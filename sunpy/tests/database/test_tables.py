from collections import Hashable
from datetime import datetime
import os.path

from sunpy.database.tables import FitsHeaderEntry, Tag, DatabaseEntry
from sunpy.net.vso import VSOClient
import sunpy

import pytest


def test_fits_header_entry_equality():
    assert FitsHeaderEntry('key', 'value') == FitsHeaderEntry('key', 'value')
    assert not (FitsHeaderEntry('key', 'value') == FitsHeaderEntry('k', 'v'))


def test_fits_header_entry_inequality():
    assert FitsHeaderEntry('key', 'value') != FitsHeaderEntry('k', 'v')
    assert not (FitsHeaderEntry('k', 'v') != FitsHeaderEntry('k', 'v'))


def test_fits_header_entry_hashability():
    assert not isinstance(FitsHeaderEntry('key', 'value'), Hashable)


def test_tag_equality():
    assert Tag('abc') == Tag('abc')
    assert not (Tag('abc') == Tag('xyz'))


def test_tag_inequality():
    assert Tag('abc') != Tag('xyz')
    assert not (Tag('abc') != Tag('abc'))


def test_tag_hashability():
    assert not isinstance(Tag(''), Hashable)


@pytest.mark.slow
def test_entry_from_qr_block():
    client = VSOClient()
    qr = client.query_legacy('2001/1/1', '2001/1/2', instrument='EIT')
    entry = DatabaseEntry.from_query_result_block(qr[0])
    expected_entry = DatabaseEntry(
        source='SOHO', provider='SDAC', physobs='intensity',
        fileid='/archive/soho/private/data/processed/eit/lz/2001/01/efz20010101.010014',
        observation_time_start=datetime(2001, 1, 1, 1, 0, 14),
        observation_time_end=datetime(2001, 1, 1, 1, 0, 21),
        instrument='EIT', size=2059.0, waveunit='Angstrom', wavemin=171.0,
        wavemax=171.0)
    assert entry == expected_entry


def test_add_fits_header_entries_from_file():
    entry = DatabaseEntry()
    assert entry.fits_header_entries == []
    entry.add_fits_header_entries_from_file(sunpy.RHESSI_EVENT_LIST)
    assert entry.fits_header_entries == [
        FitsHeaderEntry('SIMPLE', True),
        FitsHeaderEntry('BITPIX', 8),
        FitsHeaderEntry('NAXIS', 0),
        FitsHeaderEntry('EXTEND', True),
        FitsHeaderEntry('DATE', '2011-09-13T15:37:38'),
        FitsHeaderEntry('ORIGIN', 'RHESSI'),
        FitsHeaderEntry('OBSERVER', 'Unknown'),
        FitsHeaderEntry('TELESCOP', 'RHESSI'),
        FitsHeaderEntry('INSTRUME', 'RHESSI'),
        FitsHeaderEntry('OBJECT', 'Sun'),
        FitsHeaderEntry('DATE_OBS', '2002-02-20T11:06:00.000'),
        FitsHeaderEntry('DATE_END', '2002-02-20T11:06:43.330'),
        FitsHeaderEntry('TIME_UNI', 1),
        FitsHeaderEntry('ENERGY_L', 25.0),
        FitsHeaderEntry('ENERGY_H', 40.0),
        FitsHeaderEntry('TIMESYS', '1979-01-01T00:00:00'),
        FitsHeaderEntry('TIMEUNIT', 'd')]


def test_add_tags_no_params():
    with pytest.raises(TypeError):
        DatabaseEntry().add_tags()


def test_add_tags():
    entry = DatabaseEntry()
    assert entry.tags == []
    entry.add_tags('foo', 'bar', 'baz')
    assert entry.tags == [Tag('foo'), Tag('bar'), Tag('baz')]
