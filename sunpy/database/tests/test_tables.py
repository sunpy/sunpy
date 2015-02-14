# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

from collections import Hashable
from datetime import datetime

import pytest

from astropy import units as u

from sunpy.database import Database
from sunpy.database.tables import FitsHeaderEntry, FitsKeyComment, Tag,\
    DatabaseEntry, entries_from_query_result, entries_from_dir,\
    entries_from_file, display_entries, WaveunitNotFoundError
from sunpy.net import vso
from sunpy.data.test import rootdir as testdir
from sunpy.data.test.waveunit import waveunitdir, MQ_IMAGE
from sunpy.tests.helpers import skip_windows
import pytest
import os

RHESSI_IMAGE = os.path.join(testdir, 'hsi_image_20101016_191218.fits')
EIT_195_IMAGE = os.path.join(testdir, 'EIT/efz20040301.000010_s.fits')

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
        vso.attrs.Provider('SSC'), vso.attrs.Wave(10 * u.kHz, 160 * u.kHz))[0]


@pytest.fixture
def qr_block_with_kev_unit():
    return vso.VSOClient().query(
        vso.attrs.Time((2011, 9, 20, 1), (2011, 9, 20, 2)),
        vso.attrs.Instrument('RHESSI'))[0]


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


@pytest.mark.online
def test_entry_from_qr_block(query_result):
    entry = DatabaseEntry._from_query_result_block(query_result[0])
    expected_entry = DatabaseEntry(
        source='SOHO', provider='SDAC', physobs='intensity',
        fileid='/archive/soho/private/data/processed/eit/lz/2001/01/efz20010101.010014',
        observation_time_start=datetime(2001, 1, 1, 1, 0, 14),
        observation_time_end=datetime(2001, 1, 1, 1, 0, 21),
        instrument='EIT', size=2059.0, wavemin=17.1, wavemax=17.1)
    assert entry == expected_entry


@pytest.mark.online
def test_entry_from_qr_block_with_missing_physobs(qr_block_with_missing_physobs):
    entry = DatabaseEntry._from_query_result_block(qr_block_with_missing_physobs)
    expected_entry = DatabaseEntry(
        source='STEREO_A', provider='SSC',
        fileid='swaves/2013/swaves_average_20130805_a_hfr.dat',
        observation_time_start=datetime(2013, 8, 5),
        observation_time_end=datetime(2013, 8, 6), instrument='SWAVES',
        size=3601.08, wavemin=2398339664000.0, wavemax=18737028625.0)
    assert entry == expected_entry


@pytest.mark.online
def test_entry_from_qr_block_kev(qr_block_with_kev_unit):
    # See issue #766.
    entry = DatabaseEntry._from_query_result_block(qr_block_with_kev_unit)
    assert entry.source == 'RHESSI'
    assert entry.provider == 'LSSP'
    assert entry.fileid == '/hessidata/2011/09/20/hsi_20110920_010920'
    assert entry.observation_time_start == datetime(2011, 9, 20, 1, 9, 20)
    assert entry.observation_time_end == datetime(2011, 9, 20, 2, 27, 40)
    assert entry.instrument == 'RHESSI'
    assert entry.size == -1
    assert round(entry.wavemin, 3) == 0.413
    assert round(entry.wavemax, 7) == 0.0000729


def test_entries_from_file():
    entries = list(entries_from_file(MQ_IMAGE))
    assert len(entries) == 1
    entry = entries[0]
    assert len(entry.fits_header_entries) == 31
    expected_fits_header_entries = [
        FitsHeaderEntry('SIMPLE', True),
        FitsHeaderEntry('BITPIX', 16),
        FitsHeaderEntry('NAXIS', 2),
        FitsHeaderEntry('NAXIS1', 1500),
        FitsHeaderEntry('NAXIS2', 1340),
        FitsHeaderEntry('CONTACT', 'Isabelle.Buale@obspm.fr'),
        FitsHeaderEntry('DATE_OBS', '2013-08-12T08:42:53.000'),
        FitsHeaderEntry('DATE_END', '2013-08-12T08:42:53.000'),
        FitsHeaderEntry('FILENAME', 'mq130812.084253.fits'),
        FitsHeaderEntry('INSTITUT', 'Observatoire de Paris'),
        FitsHeaderEntry('INSTRUME', 'Spectroheliograph'),
        FitsHeaderEntry('OBJECT', 'FS'),
        FitsHeaderEntry('OBS_MODE', 'SCAN'),
        FitsHeaderEntry('PHYSPARA', 'Intensity'),
        FitsHeaderEntry('NBREG', 1),
        FitsHeaderEntry('NBLAMBD', 1),
        FitsHeaderEntry('WAVELNTH', 6563),
        FitsHeaderEntry('WAVEUNIT', 'angstrom'),
        FitsHeaderEntry('POLARANG', 0),
        FitsHeaderEntry('THEMISFF', 3),
        FitsHeaderEntry('LONGTRC', 258.78),
        FitsHeaderEntry('LONGCARR', 258.78),
        FitsHeaderEntry('LONGITUD', 258.78),
        FitsHeaderEntry('LATITUD', 6.50107),
        FitsHeaderEntry('LATIRC', 6.50107),
        FitsHeaderEntry('INDLAMD', 1),
        FitsHeaderEntry('INDREG', 1),
        FitsHeaderEntry('SEQ_IND', 1),
        FitsHeaderEntry('SVECTOR', 0),
        FitsHeaderEntry('COMMENT', ''),
        FitsHeaderEntry('HISTORY', '')]
    assert entry.fits_header_entries == expected_fits_header_entries
    assert entry.fits_key_comments == [
        FitsKeyComment('SIMPLE', 'Written by IDL:  Mon Aug 12 08:48:08 2013'),
        FitsKeyComment('BITPIX', 'Integer*2 (short integer)')]
    assert entry.instrument == 'Spectroheliograph'
    assert entry.observation_time_start == datetime(2013, 8, 12, 8, 42, 53)
    assert entry.observation_time_end == datetime(2013, 8, 12, 8, 42, 53)
    assert round(entry.wavemin, 1) == 656.3
    assert round(entry.wavemax, 1) == 656.3
    assert entry.path == MQ_IMAGE


def test_entries_from_file_withoutwaveunit():
    # does not raise `WaveunitNotFoundError`, because no wavelength information
    # is present in this file
    entries_from_file(RHESSI_IMAGE).next()
    with pytest.raises(WaveunitNotFoundError):
        entries_from_file(EIT_195_IMAGE).next()


def test_entries_from_dir():
    entries = list(entries_from_dir(waveunitdir))
    assert len(entries) == 4
    for entry, filename in entries:
        if filename.endswith('na120701.091058.fits'):
            break
    assert entry.path in (os.path.join(waveunitdir, filename), filename)
    assert filename.startswith(waveunitdir)
    assert len(entry.fits_header_entries) == 42
    assert entry.fits_header_entries == [
        FitsHeaderEntry('SIMPLE', True),
        FitsHeaderEntry('BITPIX', -32),
        FitsHeaderEntry('NAXIS', 3),
        FitsHeaderEntry('NAXIS1', 256),
        FitsHeaderEntry('NAXIS2', 256),
        FitsHeaderEntry('NAXIS3', 1),
        FitsHeaderEntry('DATE', '27-OCT-82'),
        FitsHeaderEntry('DATE-OBS', '2012-07-01'),
        FitsHeaderEntry('DATE_OBS', '2012-07-01T09:10:58.200Z'),
        FitsHeaderEntry('DATE_END', '2012-07-01T09:10:58.200Z'),
        FitsHeaderEntry('WAVELNTH', 1.98669),
        FitsHeaderEntry('WAVEUNIT', 'm'),
        FitsHeaderEntry('PHYSPARA', 'STOKESI'),
        FitsHeaderEntry('OBJECT', 'FS'),
        FitsHeaderEntry('OBS_TYPE', 'RADIO'),
        FitsHeaderEntry('OBS_MODE', 'IMAGE'),
        FitsHeaderEntry('LONGITUD', 0.0),
        FitsHeaderEntry('LATITUDE', 0.0),
        FitsHeaderEntry('INSTITUT', 'MEUDON'),
        FitsHeaderEntry('CMP_NAME', 'ROUTINE'),
        FitsHeaderEntry('CONTACT', ' A. KERDRAON'),
        FitsHeaderEntry('TELESCOP', 'NRH'),
        FitsHeaderEntry('INSTRUME', 'NRH2'),
        FitsHeaderEntry('FILENAME', 'nrh2_1509_h80_20120701_091058c02_i.fts'),
        FitsHeaderEntry('NRH_DATA', '2DB'),
        FitsHeaderEntry('ORIGIN', 'wrfits'),
        FitsHeaderEntry('FREQ', 150.9),
        FitsHeaderEntry('FREQUNIT', 6),
        FitsHeaderEntry('BSCALE', 1.0),
        FitsHeaderEntry('BZERO', 0.0),
        FitsHeaderEntry('BUNIT', 'K'),
        FitsHeaderEntry('EXPTIME', 1168576512),
        FitsHeaderEntry('CTYPE1', 'Solar-X'),
        FitsHeaderEntry('CTYPE2', 'Solar-Y'),
        FitsHeaderEntry('CTYPE3', 'StokesI'),
        FitsHeaderEntry('CRPIX1', 128),
        FitsHeaderEntry('CRPIX2', 128),
        FitsHeaderEntry('CDELT1', 0.015625),
        FitsHeaderEntry('CDELT2', 0.015625),
        FitsHeaderEntry('SOLAR_R', 64.0),
        FitsHeaderEntry('COMMENT', ''),
        FitsHeaderEntry('HISTORY', '')]
    assert entry.fits_key_comments == [
        FitsKeyComment('WAVEUNIT', 'in meters'),
        FitsKeyComment('NAXIS2', 'number of rows'),
        FitsKeyComment('CDELT2', 'pixel scale y, in solar radius/pixel'),
        FitsKeyComment('CRPIX1', 'SUN CENTER X, pixels'),
        FitsKeyComment('CRPIX2', 'SUN CENTER Y, pixels'),
        FitsKeyComment('SOLAR_R', 'SOLAR RADIUS, pixels'),
        FitsKeyComment('NAXIS1', 'number of columns'),
        FitsKeyComment('CDELT1', 'pixel scale x, in solar radius/pixel'),
        FitsKeyComment('NAXIS3', 'StokesI'),
        FitsKeyComment('TELESCOP', 'Nancay Radioheliograph'),
        FitsKeyComment('INSTRUME', 'Nancay 2D-images Radioheliograph'),
        FitsKeyComment('BUNIT', 'Brightness temperature'),
        FitsKeyComment('BITPIX', 'IEEE 32-bit floating point values'),
        FitsKeyComment('DATE', 'Date of file creation'),
        FitsKeyComment('FREQUNIT', 'in MHz'),
        FitsKeyComment('EXPTIME', 'in seconds')]


def test_entries_from_dir_recursively_true():
    entries = list(
        entries_from_dir(testdir, True, default_waveunit='angstrom'))
    assert len(entries) == 59
    # Older val = 31.

def test_entries_from_dir_recursively_false():
    entries = list(
        entries_from_dir(testdir, False, default_waveunit='angstrom'))
    assert len(entries) == 38


@pytest.mark.online
def test_entries_from_query_result(query_result):
    entries = list(entries_from_query_result(query_result))
    assert len(entries) == 122
    snd_entry = entries[1]
    expected_entry = DatabaseEntry(
        source='SOHO', provider='SDAC', physobs='intensity',
        fileid='/archive/soho/private/data/processed/eit/lz/2001/01/efz20010101.070014',
        observation_time_start=datetime(2001, 1, 1, 7, 0, 14),
        observation_time_end=datetime(2001, 1, 1, 7, 0, 21),
        instrument='EIT', size=2059.0, wavemin=17.1, wavemax=17.1)
    assert snd_entry == expected_entry


@pytest.mark.online
def test_entry_from_query_results_with_none_wave(qr_with_none_waves):
    # does not raise WaveunitNotFoundError because neither wavemin nor wavemax
    # are given
    list(entries_from_query_result(qr_with_none_waves))


@pytest.mark.online
def test_entry_from_query_results_with_none_wave_and_default_unit(
        qr_with_none_waves):
    entries = list(entries_from_query_result(qr_with_none_waves, 'nm'))
    assert len(entries) == 4
    assert entries == [
        DatabaseEntry(
            source='SOHO', provider='SDAC', physobs='intensity',
            fileid='/archive/soho/private/data/processed/virgo/level1/1212/HK/121222_1.H01',
            observation_time_start=datetime(2012, 12, 23, 23, 59, 3),
            observation_time_end=datetime(2012, 12, 24, 23, 59, 2),
            instrument='VIRGO', size=155.0, wavemin=None,
            wavemax=None),
        DatabaseEntry(
            source='SOHO', provider='SDAC', physobs='intensity',
            fileid='/archive/soho/private/data/processed/virgo/level1/1212/LOI/121224_1.L01',
            observation_time_end=datetime(2012, 12, 24, 23, 59, 2),
            observation_time_start=datetime(2012, 12, 23, 23, 59, 3),
            instrument='VIRGO', size=329.0, wavemin=None,
            wavemax=None),
        DatabaseEntry(
            source='SOHO', provider='SDAC', physobs ='intensity',
            fileid='/archive/soho/private/data/processed/virgo/level1/1212/SPM/121222_1.S02',
            observation_time_start=datetime(2012, 12, 23, 23, 59, 3),
            observation_time_end=datetime(2012, 12, 24, 23, 59, 2),
            instrument='VIRGO', size=87.0, wavemin=None,
            wavemax=None),
        DatabaseEntry(
            source='SOHO', provider='SDAC', physobs='intensity',
            fileid='/archive/soho/private/data/processed/virgo/level1/1212/DIARAD/121222_1.D01',
            observation_time_start=datetime(2012, 12, 24, 0, 1, 58),
            observation_time_end=datetime(2012, 12, 25, 0, 1, 57),
            instrument='VIRGO', size=14.0, wavemin=None,
            wavemax=None)]


def test_display_entries_missing_entries():
    with pytest.raises(TypeError):
        display_entries([], ['some', 'columns'])


def test_display_entries_empty_db():
    with pytest.raises(TypeError):
        display_entries(Database('sqlite:///'), ['id'])


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
            instrument='EIT', size=259.0, wavemin=171.0,
            wavemax=171.0, tags=[Tag('foo'), Tag('bar')]),
        DatabaseEntry(
            id=2, source='GONG', provider='NSO', physobs='LOS_velocity',
            fileid='pptid=11010...',
            observation_time_start=datetime(2010, 1, 1, 0, 59),
            observation_time_end=datetime(2010, 1, 1, 1),
            download_time=datetime(2014, 6, 15, 3, 42, 55, 123456),
            instrument='Merged gong', size=944.0,
            wavemin=6768.0, wavemax=6768.0, starred=True)]
    columns = [
        'id', 'source', 'provider', 'physobs', 'fileid', 'download_time',
        'observation_time_start', 'instrument', 'size',
        'wavemin', 'path', 'starred', 'tags']
    table = display_entries(entries, columns)
    filedir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(filedir,'test_table.txt'), 'r') as f:
        stored_table = f.read()
    assert table.strip() == stored_table.strip()
