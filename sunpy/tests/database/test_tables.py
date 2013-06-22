from collections import Hashable
from datetime import datetime
import os.path

from sunpy.database.tables import FitsHeaderEntry, Tag, DatabaseEntry
from sunpy.net.vso import VSOClient

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
    fitsfile = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            'lyra_20110810-000000_lev2_std.fits'))
    entry.add_fits_header_entries_from_file(fitsfile)
    assert entry.fits_header_entries == [
        FitsHeaderEntry('SIMPLE', True),
        FitsHeaderEntry('BITPIX', 8),
        FitsHeaderEntry('NAXIS', 0),
        FitsHeaderEntry('EXTEND', True),
        FitsHeaderEntry('ORIGIN', 'ROB'),
        FitsHeaderEntry('TELESCOP', 'PROBA2'),
        FitsHeaderEntry('INSTRUME', 'LYRA'),
        FitsHeaderEntry('OBJECT', 'EUV solar irrad'),
        FitsHeaderEntry('OBS_MODE', 'standard'),
        FitsHeaderEntry('DATE', '2012-10-12'),
        FitsHeaderEntry('DATE-OBS', '2011-08-10T00:00:00.020000'),
        FitsHeaderEntry('DATE-END', '2011-08-10T23:59:59.982999'),
        FitsHeaderEntry('DATASRC', 'Redu'),
        FitsHeaderEntry('LEVEL', '2'),
        FitsHeaderEntry('ALGOR_V', 'EDG=2.1  BSDG=0.8'),
        FitsHeaderEntry('FILENAME', 'lyra_20110810-000000_lev2_std.fits')]


def test_add_tags_no_params():
    with pytest.raises(TypeError):
        DatabaseEntry().add_tags()


def test_add_tags():
    entry = DatabaseEntry()
    assert entry.tags == []
    entry.add_tags('foo', 'bar', 'baz')
    assert entry.tags == [Tag('foo'), Tag('bar'), Tag('baz')]
