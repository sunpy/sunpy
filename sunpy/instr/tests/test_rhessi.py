# -*- coding: utf-8 -*-
"""
Unit tests for `sunpy.instr.rhessi`
"""
import os
import socket
import warnings
from datetime import datetime
import textwrap

import sunpy.map
import sunpy.data.test
import sunpy.instr.rhessi as rhessi
from sunpy.extern.six.moves.urllib.error import URLError

import numpy as np
import pytest
import mock

testpath = sunpy.data.test.rootdir


def test_backprojection():
    amap = rhessi.backprojection(os.path.join(testpath, 'hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits'))
    assert isinstance(amap, sunpy.map.GenericMap)
    assert amap.date == datetime(2002, 2, 20, 11, 6, 21)


def test_get_obssumm_dbase_file():
    with pytest.raises(ValueError):
        rhessi.get_obssumm_dbase_file(['2002/01/01', '2002/04/01'])


@pytest.mark.remote_data
def test_get_obssum_filename():
    file_name = rhessi.get_obssum_filename(('2011/04/04', '2011/04/05'))
    # Irregardless of mirror server the osbsum file name should match
    assert file_name[0].split('metadata/catalog/')[1] == 'hsi_obssumm_20110404_042.fits'


@pytest.mark.remote_data
def test_parse_obssum_dbase_file():
    file = rhessi.get_obssumm_dbase_file(('2011/04/04', '2011/04/05'))
    obssum = rhessi.parse_obssumm_dbase_file(file[0])

    assert obssum['filename'][0] == 'hsi_obssumm_20110401_043.fit'
    assert obssum['filename'][-1] == 'hsi_obssumm_20110430_029.fit'

    assert obssum['orb_st'][0] == 0
    assert obssum['orb_st'][-1] == 0

    assert obssum['orb_end'][0] == 0
    assert obssum['orb_end'][-1] == 0

    assert obssum['start_time'][0] == datetime(2011, 4, 1, 0, 0, 0)
    assert obssum['start_time'][-1] == datetime(2011, 4, 30, 0, 0, 0)

    assert obssum['end_time'][0] == datetime(2011, 4, 2, 0, 0, 0)
    assert obssum['end_time'][-1] == datetime(2011, 5, 1, 0, 0, 0)

    assert obssum['status_flag'][0] == 0
    assert obssum['status_flag'][-1] == 0

    assert obssum['npackets'][0] == 0
    assert obssum['npackets'][-1] == 0


@pytest.mark.remote_data
def test_get_parse_obssum_file():
    f = rhessi.get_obssumm_file(('2011/04/04', '2011/04/05'))
    header, _data = rhessi.parse_obssumm_file(f[0])
    assert header.get('DATE_OBS') == '2011-04-04T00:00:00.000'
    assert header.get('DATE_END') == '2011-04-05T00:00:00.000'
    assert header.get('TELESCOP') == 'HESSI'


def test_uncompress_countrate():
    # Should only accept bytearr (uncompressed counts must be 0 - 255)
    with pytest.raises(ValueError):
        rhessi.uncompress_countrate(np.array([-1, 300]))

    counts = rhessi.uncompress_countrate(np.array([0, 128, 255]))

    # Valid min, max
    assert counts[0] == 0
    assert counts[2] == 1015792

    # Random test value
    assert counts[1] == 4080


@pytest.fixture
def unsupported_start_time():
    """
    RHESSI summary files are not available for before 2002-02-01
    """
    return sunpy.time.TimeRange(("2002/01/31", "2002/02/03"))


@pytest.fixture
def one_day_timerange():
    return sunpy.time.TimeRange(("2016/01/15", "2016/01/16"))


@pytest.fixture
def two_days_timerange():
    return sunpy.time.TimeRange(("2016/01/15", "2016/01/17"))


# Test `rhessi.get_base_url()`

@mock.patch('sunpy.instr.rhessi.urlopen', return_value=None)
def test_get_base_url(mock_urlopen):
    """
    Success case, can successfully 'ping' first data_server
    """
    assert rhessi.get_base_url() == rhessi.data_servers[0]


@mock.patch('sunpy.instr.rhessi.urlopen', side_effect=URLError(''))
def test_get_base_url_on_urlerror(mock_urlopen):
    """
    If all tested URLs raise `URLError`, then raise an `IOError`
    """
    with pytest.raises(IOError):
        rhessi.get_base_url()


@mock.patch('sunpy.instr.rhessi.urlopen', side_effect=socket.timeout)
def test_get_base_url_on_timeout(mock_urlopen):
    """
    If all tested data servers timeout, then raise an `IOError`
    """
    with pytest.raises(IOError):
        rhessi.get_base_url()


# Test `rhessi.get_obssumm_dbase_file(...)`

def test_get_obssumm_dbase_file_with_unsupported_start_time(unsupported_start_time):
    """
    RHESSI summary files are not available for before 2002-02-01, ensure
    `ValueError` is raised.
    """
    with pytest.raises(ValueError):
        rhessi.get_obssumm_dbase_file(unsupported_start_time)


@mock.patch('sunpy.instr.rhessi.urlretrieve', return_value=None)
@mock.patch('sunpy.instr.rhessi.get_base_url', return_value='http://www.example.com')
def test_get_obssumm_dbase_file_timerange_warnings(mock_get_base_url,
                                                   mock_urlretrieve,
                                                   one_day_timerange,
                                                   two_days_timerange):
    """
    With only one day in the time range, the `UserWarning` should *not* be issued.
    With two days in the time range, the `UserWarning` should be issued.

    Only interested in the presence/absence of the `UserWarning`, *not* the return
    value of `get_obssumm_dbase_file`.
    """
    with warnings.catch_warnings(record=True) as caught_warnings:
        rhessi.get_obssumm_dbase_file(one_day_timerange)
        assert len(caught_warnings) == 0

    with pytest.warns(UserWarning):
        rhessi.get_obssumm_dbase_file(two_days_timerange)


@mock.patch('sunpy.instr.rhessi.urlretrieve', return_value=None)
@mock.patch('sunpy.instr.rhessi.get_base_url', return_value='http://www.example.com')
def test_get_obssumm_dbase_file_build_correct_url(mock_get_base_url, mock_urlretrieve,
                                                  one_day_timerange):
    """
    This test ensures that we build the correct url which is then used
    to get the database file.
    """
    rhessi.get_obssumm_dbase_file(one_day_timerange)
    mock_urlretrieve.assert_called_with(
        'http://www.example.com/dbase/hsi_obssumm_filedb_201601.txt')


# Test `rhessi.parse_obssumm_dbase_file(...)`

def hessi_data():
    return textwrap.dedent("""\
         HESSI Filedb File:
         Created: 1972-04-14T12:41:26.000
         Number of Files:           2
                             Filename  Orb_st Orb_end         Start_time           End_time Status_flag    Npackets Drift_start   Drift_end Data source
         hsi_obssumm_19721101_139.fit       7       8 01-Nov-72 00:00:00 02-Nov-72 00:00:00           3           2       0.000       0.000
         hsi_obssumm_19721102_144.fit       9      10 02-Nov-72 00:00:00 03-Nov-72 00:00:00           4           1       0.000       0.000
         """).splitlines()


def test_parse_obssumm_dbase_file():
    """
    Ensure that all required data are extracted from the RHESSI
    observing summary database file mocked in `hessi_data()`
    """
    mock_file = mock.mock_open()
    mock_file.return_value.__iter__.return_value = hessi_data()

    dbase_data = {}
    with mock.patch('sunpy.instr.rhessi.open', mock_file, create=True):
        dbase_data = rhessi.parse_obssumm_dbase_file(None)

    assert len(dbase_data.keys()) == 7

    # verify each of the 7 fields
    assert dbase_data['filename'] == ['hsi_obssumm_19721101_139.fit',
                                      'hsi_obssumm_19721102_144.fit']
    assert dbase_data['orb_st'] == [7, 9]
    assert dbase_data['orb_end'] == [8, 10]
    assert dbase_data['start_time'] == [datetime(1972, 11, 1, 0, 0), datetime(1972, 11, 2, 0, 0)]
    assert dbase_data['end_time'] == [datetime(1972, 11, 2, 0, 0), datetime(1972, 11, 3, 0, 0)]
    assert dbase_data['status_flag'] == [3, 4]
    assert dbase_data['npackets'] == [2, 1]


# Test `rhessi.get_obssum_filename(...)`

def parsed_dbase():
    """
    The result of calling `parse_obssumm_dbase_file(...)` on
    https://hesperia.gsfc.nasa.gov/hessidata/dbase/hsi_obssumm_filedb_200311.txt but
    only using the first two rows of data.
    """
    return {'filename': ['hsi_obssumm_20031101_139.fit', 'hsi_obssumm_20031102_144.fit'],
            'orb_st': [0, 0],
            'orb_end': [0, 0],
            'start_time': [datetime(2003, 11, 1, 0, 0), datetime(2003, 11, 2, 0, 0)],
            'end_time': [datetime(2003, 11, 2, 0, 0), datetime(2003, 11, 3, 0, 0)],
            'status_flag': [0, 0],
            'npackets': [0, 0]}


@mock.patch('sunpy.instr.rhessi.get_base_url', return_value='http://www.example.com')
@mock.patch('sunpy.instr.rhessi.parse_obssumm_dbase_file', return_value=parsed_dbase())
@mock.patch('sunpy.instr.rhessi.get_obssumm_dbase_file', return_value=('', {}))
def test_get_obssum_filename_one_day(mock_get_obssumm_dbase_file,
                                     mock_parse_obssumm_dbase_file,
                                     mock_get_base_url):
    """
    Given a time range of one day, make sure we get one days data back, i.e. one file.
    """
    filename = rhessi.get_obssum_filename(('2003-11-01', '2003-11-02'))

    assert len(filename) == 1
    assert filename[0] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031101_139.fits'


@mock.patch('sunpy.instr.rhessi.get_base_url', return_value='http://www.example.com')
@mock.patch('sunpy.instr.rhessi.parse_obssumm_dbase_file', return_value=parsed_dbase())
@mock.patch('sunpy.instr.rhessi.get_obssumm_dbase_file', return_value=('', {}))
def test_get_obssum_filename_two_days(mock_get_obssumm_dbase_file,
                                      mock_parse_obssumm_dbase_file,
                                      mock_get_base_url):
    """
    Given a time range of two days, make sure we get two files back, one
    for each day.
    """
    filenames = rhessi.get_obssum_filename(('2003-11-01', '2003-11-03'))

    assert len(filenames) == 2
    assert filenames[0] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031101_139.fits'
    assert filenames[1] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031102_144.fits'


# Test `rhessi.get_obssumm_file(...)`

@mock.patch('sunpy.instr.rhessi.urlretrieve', return_value=None)
@mock.patch('sunpy.instr.rhessi.get_obssum_filename', return_value=['dummy.fits'])
def test_get_obssumm_file_timerange_warnings(mock_get_obssum_filename,
                                             mock_urlretrieve,
                                             one_day_timerange,
                                             two_days_timerange):
    """
    With only one day in the time range, the `UserWarning` should *not* be issued.
    With two days in the time range, the `UserWarning` should be issued.

    Only interested in the presence/absence of the `UserWarning`, *not* the return
    value of `get_obssumm_file`.
    """
    with warnings.catch_warnings(record=True) as caught_warnings:
        rhessi.get_obssumm_file(one_day_timerange)
        assert len(caught_warnings) == 0

    with pytest.warns(UserWarning):
        rhessi.get_obssumm_file(two_days_timerange)


@mock.patch('sunpy.instr.rhessi.urlretrieve', return_value=None)
@mock.patch('sunpy.instr.rhessi.get_obssum_filename', return_value=['http://first.org/file.txt',
                                                                    'http://second.com/f.txt'])
def test_get_obssumm_file_timerange_warnings(mock_get_obssum_filename,
                                             mock_urlretrieve,
                                             two_days_timerange):
    """
    Even though `get_obssum_filename` returns more than one file. Make sure that
    `get_obssumm_file` only attempts to retrieve the first one.
    """
    rhessi.get_obssumm_file(two_days_timerange)

    mock_urlretrieve.assert_called_with('http://first.org/file.txt')


# Test `rhessi._build_energy_bands(...)`

@pytest.fixture
def raw_bands():
    return ['3 - 6', '6 - 12', '12 - 25', '25 - 50', '50 - 100', '100 - 300', '300 - 800',
            '800 - 7000', '7000 - 20000']


def test__build_energy_bands_no_match(raw_bands):
    """
    If an energy unit cannot be found in the `label` then raise
    a `ValueError`
    """
    with pytest.raises(ValueError):
        rhessi._build_energy_bands(label='Energy bands GHz', bands=raw_bands)


def test__build_energy_bands(raw_bands):
    """
    Success case.
    """
    built_ranges = rhessi._build_energy_bands(label='Energy bands (keV)', bands=raw_bands)

    assert built_ranges == ['3 - 6 keV', '6 - 12 keV', '12 - 25 keV', '25 - 50 keV',
                            '50 - 100 keV', '100 - 300 keV', '300 - 800 keV',
                            '800 - 7000 keV', '7000 - 20000 keV']


# Test `rhessi._check_one_day(...)`

def test__check_one_day_with_warning(two_days_timerange):
    """
    Issue a warning if TimeRange > one day
    """
    with pytest.warns(UserWarning):
        rhessi._check_one_day(two_days_timerange)


def test__check_one_day(one_day_timerange):
    """
    Success case, time range is one day - no warnings should be raised
    """

    with warnings.catch_warnings(record=True) as caught_warnings:
        rhessi._check_one_day(one_day_timerange)
        assert len(caught_warnings) == 0
