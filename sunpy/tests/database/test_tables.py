from collections import Hashable
from datetime import datetime
import os

from sunpy.database.tables import FitsHeaderEntry, Tag, DatabaseEntry,\
    entries_from_query_result, entries_from_path
from sunpy.net.vso import VSOClient
from sunpy.data.test import rootdir as testdir
import sunpy

import pytest


@pytest.fixture
def query_result():
    client = VSOClient()
    return client.query_legacy('2001/1/1', '2001/1/2', instrument='EIT')


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
    assert isinstance(Tag(''), Hashable)


@pytest.mark.slow
def test_entry_from_qr_block(query_result):
    entry = DatabaseEntry.from_query_result_block(query_result[0])
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
    assert entry.instrument == 'RHESSI'
    assert entry.observation_time_start == datetime(2002, 02, 20, 11, 6, 0, 0)


def test_add_fits_header_entries_from_file_wavelength():
    entry = DatabaseEntry()
    assert entry.fits_header_entries == []
    path = os.path.join(testdir, 'EIT', 'efz20040301.020010_s.fits')
    entry.add_fits_header_entries_from_file(path)
    assert FitsHeaderEntry('WAVELNTH', 195) in entry.fits_header_entries
    assert entry.wavemin == 195.0
    assert entry.wavemax == 195.0


def test_entries_from_path():
    entries = list(entries_from_path(os.path.join(testdir, 'EIT')))
    assert len(entries) == 13
    assert entries[0].fits_header_entries == [
        FitsHeaderEntry('SIMPLE', True),
        FitsHeaderEntry('BITPIX', -64),
        FitsHeaderEntry('NAXIS', 2),
        FitsHeaderEntry('NAXIS1', 128),
        FitsHeaderEntry('NAXIS2', 128),
        FitsHeaderEntry('DATE_OBS', '2004-03-01T02:00:10.642Z'),
        FitsHeaderEntry('CMP_NO', 1),
        FitsHeaderEntry('DSUN_OBS', 18501183493.45421),
        FitsHeaderEntry('SOLAR_R', 372.27),
        FitsHeaderEntry('CDELT1', 21.04),
        FitsHeaderEntry('CDELT2', 21.04),
        FitsHeaderEntry('OBJECT', 'full FOV'),
        FitsHeaderEntry('DATE-OBS', '2004-03-01T02:00:10.642Z'),
        FitsHeaderEntry('CFTEMP', -66.71),
        FitsHeaderEntry('DATE', '2004-03-01'),
        FitsHeaderEntry('EXPMODE', 'backside'),
        FitsHeaderEntry('COMMENT', "CORRECTED DATE_OBS = '2004-03-01T01:58:31.604Z'  COMMANDED EXPOSURE TIME"),
        FitsHeaderEntry('COMMENT', " =   10.000 s  SHUTTER CLOSE TIME =    2.598 s  LINE_SYNC = 'no'  CAMERA"),
        FitsHeaderEntry('COMMENT', "_ERR = 'no'  IMAGE_OF_SEQ =                    0  READOUT PORT = 'B'  NU"),
        FitsHeaderEntry('COMMENT', 'M_LEB_PROC =                    3  LEB_PROC = 26 (no ROI)  LEB_PROC = 27'),
        FitsHeaderEntry('COMMENT', ' (no occ mask)  LEB_PROC = 12 (Rice)  BLOCKS_HORZ =   32  BLOCKS_VERT ='),
        FitsHeaderEntry('COMMENT', '  32  P1_X =           1  P2_X =        1024  P1_Y =          20  P2_Y ='),
        FitsHeaderEntry('COMMENT', '        1043  N_MISSING_BLOCKS =    0'),
        FitsHeaderEntry('COMMENT', '[\n\n\n\n                       / 284 = Fe XV, 304 = He II\n\n\n\n\n\n\n\n,   Versio'),
        FitsHeaderEntry('COMMENT', 'n 4.0, 2001 December 10]'),
        FitsHeaderEntry('CAR_ROT', 2013.0),
        FitsHeaderEntry('OBS_PROG', '195_10S_AL_1.000'),
        FitsHeaderEntry('SC_Y0', 0.0),
        FitsHeaderEntry('FILENAME', 'efz20040301.020010'),
        FitsHeaderEntry('INSTRUME', 'EIT'),
        FitsHeaderEntry('CTYPE2', 'Solar-Y'),
        FitsHeaderEntry('ORIGIN', 'Rocket Science'),
        FitsHeaderEntry('CTYPE1', 'Solar-X'),
        FitsHeaderEntry('DATASRC', 'LZ file'),
        FitsHeaderEntry('SOLAR_B0', -7.22),
        FitsHeaderEntry('CCDTEMP', 7.34),
        FitsHeaderEntry('SC_X0', 0.0),
        FitsHeaderEntry('BUNIT', 'counts / pixel'),
        FitsHeaderEntry('DETECTOR', 'EIT'),
        FitsHeaderEntry('CRVAL2', -11.70350000000008),
        FitsHeaderEntry('CRPIX1', 64.5),
        FitsHeaderEntry('CRPIX2', 64.5),
        FitsHeaderEntry('CRVAL1', 15.67479999999978),
        FitsHeaderEntry('TIME-OBS', '02:00:10'),
        FitsHeaderEntry('TELESCOP', 'SOHO'),
        FitsHeaderEntry('WAVELNTH', 195),
        FitsHeaderEntry('FILTER', 'Al +1'),
        FitsHeaderEntry('SC_ROLL', 180.0),
        FitsHeaderEntry('HEC_X', -139130208.0),
        FitsHeaderEntry('HEC_Y', 46576272.0),
        FitsHeaderEntry('HEC_Z', -94896.31),
        FitsHeaderEntry('EXPTIME', 12.598),
        FitsHeaderEntry('SCI_OBJ', 'CME WATCH 195')]


def test_entries_from_path_recursively_true():
    entries = list(entries_from_path(testdir, True))
    assert len(entries) == 15


def test_entries_from_path_recursively_false():
    entries = list(entries_from_path(testdir, False))
    assert len(entries) == 2


@pytest.mark.slow
def test_entries_from_query_result(query_result):
    entries = list(entries_from_query_result(query_result))
    assert len(entries) == 122
    snd_entry = entries[1]
    expected_entry = DatabaseEntry(
        source='SOHO', provider='SDAC', physobs='intensity',
        fileid='/archive/soho/private/data/processed/eit/lz/2001/01/efz20010101.070014',
        observation_time_start=datetime(2001, 1, 1, 7, 0, 14),
        observation_time_end=datetime(2001, 1, 1, 7, 0, 21),
        instrument='EIT', size=2059.0, waveunit='Angstrom', wavemin=171.0,
        wavemax=171.0)
    assert snd_entry == expected_entry
