# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

from collections import Hashable
from datetime import datetime
import os

from sunpy.database.tables import FitsHeaderEntry, Tag, DatabaseEntry,\
    entries_from_query_result, entries_from_path, display_entries
from sunpy.net import vso
from sunpy.data.test import rootdir as testdir
import sunpy

import pytest


@pytest.fixture
def query_result():
    client = vso.VSOClient()
    return client.query_legacy('2001/1/1', '2001/1/2', instrument='EIT')


@pytest.fixture
def qr_with_none_waves():
    return vso.VSOClient().query(
        vso.attrs.Time('20121224T120049.8', '20121224T120049.8'),
        vso.attrs.Provider('SDAC'), vso.attrs.Instrument('VIRGO'))


@pytest.fixture
def qr_block_with_missing_physobs():
    return vso.VSOClient().query(
        vso.attrs.Time('20130805T120000', '20130805T121000'),
        vso.attrs.Instrument('SWAVES'), vso.attrs.Source('STEREO_A'),
        vso.attrs.Provider('SSC'), vso.attrs.Wave(10, 160, 'kHz'))[0]


def test_fits_header_entry_equality():
    assert FitsHeaderEntry('key', 'value') == FitsHeaderEntry('key', 'value')
    assert not (FitsHeaderEntry('key', 'value') == FitsHeaderEntry('k', 'v'))


def test_fits_header_entry_inequality():
    assert FitsHeaderEntry('key', 'value') != FitsHeaderEntry('k', 'v')
    assert not (FitsHeaderEntry('k', 'v') != FitsHeaderEntry('k', 'v'))


def test_fits_header_entry_hashability():
    assert isinstance(FitsHeaderEntry('key', 'value'), Hashable)


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


@pytest.mark.slow
def test_entry_from_qr_block_with_missing_physobs(qr_block_with_missing_physobs):
    entry = DatabaseEntry.from_query_result_block(qr_block_with_missing_physobs)
    expected_entry = DatabaseEntry(
        source='STEREO_A', provider='SSC',
        fileid='swaves/2013/swaves_average_20130805_a_hfr.dat',
        observation_time_start=datetime(2013, 8, 5),
        observation_time_end=datetime(2013, 8, 6), instrument='SWAVES',
        size=3601.08, waveunit='MHz', wavemin=0.125, wavemax=16)
    assert entry == expected_entry


def test_add_fits_header_entries_from_file():
    entry = DatabaseEntry()
    assert entry.fits_header_entries == []
    entry.add_fits_header_entries_from_file(sunpy.RHESSI_EVENT_LIST)
    assert len(entry.fits_header_entries) == 20
    expected_fits_header_entries = [
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
        FitsHeaderEntry('TIMEUNIT', 'd'),
        FitsHeaderEntry('COMMENT', ''),
        FitsHeaderEntry('KEYCOMMENTS', '''  SIMPLE  Written by IDL:  Tue Sep 13 15:37:38 2011
  BITPIX  
   NAXIS  
  EXTEND  File contains extensions
    DATE  File creation date (YYYY-MM-DDThh:mm:ss UTC)
  ORIGIN  High Energy Solar Spectroscopic Imager
OBSERVER  Usually the name of the user who generated file
TELESCOP  Name of the Telescope or Mission
INSTRUME  Name of the instrument
  OBJECT  Object being observed
DATE_OBS  nominal U.T. date when integration of this
DATE_END  nominal U.T. date when integration of this
TIME_UNI  
ENERGY_L  
ENERGY_H  
 TIMESYS  Reference time in YYYY MM DD hh:mm:ss
TIMEUNIT  Unit for TIMEZERO, TSTARTI and TSTOPI'''),
        FitsHeaderEntry('HISTORY', '')]
    assert entry.fits_header_entries == expected_fits_header_entries
    assert entry.instrument == 'RHESSI'
    assert entry.observation_time_start == datetime(2002, 02, 20, 11, 6, 0, 0)
    assert entry.observation_time_end == datetime(2002, 02, 20, 11, 6, 43, 330000)


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
    first_entry, filename = entries[0]
    assert filename.startswith(os.path.join(testdir, 'EIT'))
    assert filename.endswith('.fits')
    assert len(first_entry.fits_header_entries) == 47
    assert first_entry.fits_header_entries == [
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
        FitsHeaderEntry('COMMENT', """
CORRECTED DATE_OBS = '2004-03-01T01:58:31.604Z'  COMMANDED EXPOSURE TIME =   10.000 s  SHUTTER CLOSE TIME =    2.598 s  LINE_SYNC = 'no'  CAMERA_ERR = 'no'  IMAGE_OF_SEQ =                    0  READOUT PORT = 'B'  NUM_LEB_PROC =                    3  LEB_PROC = 26 (no ROI)  LEB_PROC = 27 (no occ mask)  LEB_PROC = 12 (Rice)  BLOCKS_HORZ =   32  BLOCKS_VERT =  32  P1_X =           1  P2_X =        1024  P1_Y =          20  P2_Y =        1043  N_MISSING_BLOCKS =    0[



                       / 284 = Fe XV, 304 = He II







,   Version 4.0, 2001 December 10]""".strip()),
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
        FitsHeaderEntry('SCI_OBJ', 'CME WATCH 195'),
        FitsHeaderEntry('KEYCOMMENTS',
            '  SIMPLE  conforms to FITS standard\n'
            '  BITPIX  array data type\n'
            '   NAXIS  number of array dimensions\n'
            '  NAXIS1  \n  NAXIS2  \nDATE_OBS  \n'
            '  CMP_NO  \nDSUN_OBS  \n SOLAR_R  \n  CDELT1  \n  CDELT2  \n'
            '  OBJECT  \nDATE-OBS  \n  CFTEMP  \n    DATE  \n EXPMODE  \n'
            ' COMMENT  \n COMMENT  \n COMMENT  \n COMMENT  \n COMMENT  \n'
            ' COMMENT  \n COMMENT  \n COMMENT  \n COMMENT  \n CAR_ROT  \n'
            'OBS_PROG  \n   SC_Y0  \nFILENAME  \nINSTRUME  \n  CTYPE2  \n'
            '  ORIGIN  \n  CTYPE1  \n DATASRC  \nSOLAR_B0  \n CCDTEMP  \n'
            '   SC_X0  \n   BUNIT  \nDETECTOR  \n  CRVAL2  \n  CRPIX1  \n'
            '  CRPIX2  \n  CRVAL1  \nTIME-OBS  \nTELESCOP  \nWAVELNTH  \n'
            '  FILTER  \n SC_ROLL  \n   HEC_X  \n   HEC_Y  \n   HEC_Z  \n'
            ' EXPTIME  \n SCI_OBJ  '),
        FitsHeaderEntry('HISTORY', '')]


def test_entries_from_path_recursively_true():
    entries = list(entries_from_path(testdir, True))
    assert len(entries) == 18


def test_entries_from_path_recursively_false():
    entries = list(entries_from_path(testdir, False))
    assert len(entries) == 5


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


def test_entry_from_query_results_with_none_wave(qr_with_none_waves):
    entries = list(entries_from_query_result(qr_with_none_waves))
    assert len(entries) == 5
    assert entries == [
        DatabaseEntry(
            source='SOHO', provider='SDAC', physobs='intensity',
            fileid='/archive/soho/private/data/processed/virgo/level1/1212/HK/121222_1.H01',
            observation_time_start=datetime(2012, 12, 23, 23, 59, 3),
            observation_time_end=datetime(2012, 12, 24, 23, 59, 2),
            instrument='VIRGO', size=155.0, waveunit=None, wavemin=None,
            wavemax=None),
        DatabaseEntry(
            source='SOHO', provider='SDAC', physobs='intensity',
            fileid='/archive/soho/private/data/processed/virgo/level1/1212/LOI/121224_1.L01',
            observation_time_end=datetime(2012, 12, 24, 23, 59, 2),
            observation_time_start=datetime(2012, 12, 23, 23, 59, 3),
            instrument='VIRGO', size=329.0, waveunit=None, wavemin=None,
            wavemax=None),
        DatabaseEntry(
            source='SOHO', provider='SDAC', physobs ='intensity',
            fileid='/archive/soho/private/data/processed/virgo/level1/1212/SPM/121222_1.S02',
            observation_time_start=datetime(2012, 12, 23, 23, 59, 3),
            observation_time_end=datetime(2012, 12, 24, 23, 59, 2),
            instrument='VIRGO', size=87.0, waveunit=None, wavemin=None,
            wavemax=None),
        DatabaseEntry(
            source='SOHO', provider='SDAC', physobs='intensity',
            fileid='/archive/soho/private/data/processed/virgo/level1/1212/DIARAD/121222_1.D01',
            observation_time_start=datetime(2012, 12, 24, 0, 1, 58),
            observation_time_end=datetime(2012, 12, 25, 0, 1, 57),
            instrument='VIRGO', size=14.0, waveunit=None, wavemin=None,
            wavemax=None),
      DatabaseEntry(
            source='SOHO', provider='SDAC', physobs='intensity',
            fileid='/archive/soho/private/data/processed/virgo/level1/1212/DIARAD/121222_1.D01',
            observation_time_end=datetime(2012, 12, 25, 0, 1, 57),
            observation_time_start=datetime(2012, 12, 24, 0, 1, 58),
            instrument='VIRGO', size=14.0, waveunit=None, wavemin=None,
            wavemax=None)]


def test_display_entries_missing_entries():
    with pytest.raises(TypeError):
        display_entries([], ['some', 'columns'])


def test_display_entries_missing_columns():
    with pytest.raises(TypeError):
        display_entries([DatabaseEntry()], [])


def test_display_entries():
    entries = [
        DatabaseEntry(
            id=1, source='SOHO', provider='SDAC', physobs='intensity',
            fileid='/archive/soho/...',
            observation_time_start=datetime(2001, 1, 1, 7, 0, 14),
            observation_time_end=datetime(2001, 1, 1, 7, 0, 21),
            instrument='EIT', size=259.0, waveunit='Angstrom', wavemin=171.0,
            wavemax=171.0, tags=[Tag('foo'), Tag('bar')]),
        DatabaseEntry(
            id=2, source='GONG', provider='NSO', physobs='LOS_velocity',
            fileid='pptid=11010...',
            observation_time_start=datetime(2010, 1, 1, 0, 59),
            observation_time_end=datetime(2010, 1, 1, 1),
            instrument='Merged gong', size=944.0, waveunit='Angstrom',
            wavemin=6768.0, wavemax=6768.0, starred=True)]
    columns = [
        'id', 'source', 'provider', 'physobs', 'fileid',
        'observation_time_start', 'observation_time_end', 'instrument', 'size',
        'waveunit', 'wavemin', 'path', 'starred', 'tags']
    table = display_entries(entries, columns)
    assert table == """id source provider physobs      fileid            observation_time_start observation_time_end instrument  size  waveunit wavemin path starred tags    
-- ------ -------- -------      ------            ---------------------- -------------------- ----------  ----  -------- ------- ---- ------- ----    
1  SOHO   SDAC     intensity    /archive/soho/... 2001-01-01 07:00:14    2001-01-01 07:00:21  EIT         259.0 Angstrom 171.0   N/A  No      foo, bar
2  GONG   NSO      LOS_velocity pptid=11010...    2010-01-01 00:59:00    2010-01-01 01:00:00  Merged gong 944.0 Angstrom 6768.0  N/A  Yes     N/A     """
