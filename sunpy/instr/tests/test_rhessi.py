# -*- coding: utf-8 -*-
"""
Unit tests for `sunpy.instr.rhessi`
"""
import os
import socket
from datetime import datetime, timedelta
import numpy as np
import pytest
import mock

import sunpy.map
import sunpy.data.test
import sunpy.instr.rhessi as rhessi
from sunpy.extern.six.moves.urllib.error import URLError


TESTPATH = sunpy.data.test.rootdir
SPLIT_STR = 'metadata/catalog/'


@pytest.fixture
def before_rhessi_time():
    """
    Time before RHESSI was launched which was on 2002/02/01.
    """
    return datetime(2002, 1, 15)


@pytest.fixture
def unsupported_time_range():
    """
    RHESSI summary files are not available for before 2002-02-01.
    """
    return sunpy.time.TimeRange((before_rhessi_time - timedelta(days=2),
                                 before_rhessi_time - timedelta(days=1)))


@pytest.fixture
def one_day_timerange():
    """
    Time range which covers only one day.
    """
    return sunpy.time.TimeRange(("2016/01/15 01:00", "2016/01/15 07:00"))


@pytest.fixture
def two_days_timerange():
    """
    Time range which covers two days.
    """
    return sunpy.time.TimeRange(("2016/01/15", "2016/01/16"))


@pytest.fixture
def cross_month_timerange():
    """
    Time range which crosses a month boundary. Dbase files are monthly
    therefore this is to make sure that two dbase files are returned.
    """
    return sunpy.time.TimeRange(("2016/01/25", "2016/02/05"))


def test_backprojection():
    """
    Test that backprojection returns a map with the expected time.
    """
    test_filename = 'hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits'
    amap = rhessi.backprojection(os.path.join(TESTPATH, test_filename))
    assert isinstance(amap, sunpy.map.GenericMap)
    assert amap.date == datetime(2002, 2, 20, 11, 6, 21)


def test_get_obssumm_dbase_before_rhessi(before_rhessi_time):
    """
    Test that an error is returned if asking for a time before launch.
    """
    with pytest.raises(ValueError):
        rhessi.get_observing_summary_dbase_file(before_rhessi_time)


@pytest.mark.remote_data
def test_get_observing_summary_filename_one_day(one_day_timerange):
    """
    Test that one file is returned with the expected non-changing
    part of the filename.
    """
    file_names = rhessi.get_observing_summary_filename(one_day_timerange)
    # Irregardless of mirror server the observing summary filename should match
    assert len(file_names) == np.ceil(one_day_timerange.days.value)
    assert file_names[0].split(SPLIT_STR)[1][0:20] == 'hsi_obssumm_20110404'


@pytest.mark.remote_data
def test_get_observing_summary_filename_two_day(two_days_timerange):
    """
    Test that two file are returned with the expected non-changing
    part of the filenames.
    """
    file_names = rhessi.get_observing_summary_filename(two_days_timerange)
    # Irregardless of mirror server the obssumm file name should match
    assert len(file_names) == np.ceil(two_days_timerange.days.value)
    assert file_names[0].split(SPLIT_STR)[1][0:20] == 'hsi_obssumm_20110404'
    assert file_names[1].split(SPLIT_STR)[1][0:20] == 'hsi_obssumm_20110404'


@pytest.mark.remote_data
def test_get_observing_summary_filename_cross_month(cross_month_timerange):
    """
    Test that crossing a month returns the right number of files.
    part of the filename.
    """
    file_names = rhessi.get_observing_summary_filename(cross_month_timerange)
    # Irregardless of mirror server the obssumm file name should match
    assert len(file_names) == np.ceil(cross_month_timerange.days.value)


@pytest.mark.remote_data
def test_parse_observing_summary_dbase_file(one_day_timerange):
    """
    Test that we get the observing summary dbase file with the content
    we expect.
    """
    file = rhessi.get_observing_summary_filename(one_day_timerange)
    obssum = rhessi.parse_observing_summary_dbase_file(file[0])

    assert obssum['filename'][0][0:20] == 'hsi_obssumm_20110401'
    assert obssum['filename'][-1][0:20] == 'hsi_obssumm_20110430'

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
def test_get_parse_observing_summary_file(one_day_timerange):
    """
    Test that the contents of the header is what we expect.
    """
    f = rhessi.get_observing_summary_filename(one_day_timerange)
    header, _data = rhessi.parse_observing_summary_file(f[0])
    assert header.get('DATE_OBS') == '2011-04-04T00:00:00.000'
    assert header.get('DATE_END') == '2011-04-05T00:00:00.000'
    assert header.get('TELESCOP') == 'HESSI'


def test_uncompress_countrate():
    """
    Test that function fails if given uncompressed counts out of range.
    """
    # Should only accept bytearr (uncompressed counts must be 0 - 255)
    with pytest.raises(ValueError):
        rhessi.uncompress_countrate(np.array([-1, 300]))

    counts = rhessi.uncompress_countrate(np.array([0, 128, 255]))

    # Valid min, max
    assert counts[0] == 0
    assert counts[2] == 1015792

    # Random test value
    assert counts[1] == 4080


# Test `rhessi.get_base_url()`

@mock.patch('sunpy.instr.rhessi.urlopen', return_value=None)
def test_get_base_url(mockurlopen):
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

# Test `rhessi.parse_obssumm_dbase_file(...)`

def hessi_data():
    """Data expected in test file."""
    return """HESSI Filedb File:
Created: 1972-04-14T12:41:26.000
Number of Files:           2
                    Filename  Orb_st Orb_end         Start_time           End_time Status_flag    Npackets Drift_start   Drift_end Data source
hsi_obssumm_19721101_139.fit       7       8 01-Nov-72 00:00:00 02-Nov-72 00:00:00           3           2       0.000       0.000
hsi_obssumm_19721102_144.fit       9      10 02-Nov-72 00:00:00 03-Nov-72 00:00:00           4           1       0.000       0.000
""".splitlines()


def test_parse_observing_summary_dbase_file_mock():
    """
    Ensure that all required data are extracted from the RHESSI
    observing summary database file mocked in `hessi_data()`
    """
    mock_file = mock.mock_open()
    mock_file.return_value.__iter__.return_value = hessi_data()

    dbase_data = {}
    with mock.patch('sunpy.instr.rhessi.open', mock_file, create=True):
        dbase_data = rhessi.parse_observing_summary_dbase_file(None)

    assert len(dbase_data.keys()) == 7

    # verify each of the 7 fields
    assert dbase_data['filename'] == ['hsi_obssumm_19721101_139.fit',
                                      'hsi_obssumm_19721102_144.fit']
    assert dbase_data['orb_st'] == [7, 9]
    assert dbase_data['orb_end'] == [8, 10]
    assert dbase_data['start_time'] == [datetime(1972, 11, 1, 0, 0),
                                        datetime(1972, 11, 2, 0, 0)]
    assert dbase_data['end_time'] == [datetime(1972, 11, 2, 0, 0),
                                      datetime(1972, 11, 3, 0, 0)]
    assert dbase_data['status_flag'] == [3, 4]
    assert dbase_data['npackets'] == [2, 1]


# Test `rhessi.get_observing_summary_filename(...)`

def parsed_dbase():
    """
    The result of calling `parse_obssumm_dbase_file(...)` on
    https://hesperia.gsfc.nasa.gov/hessidata/dbase/hsi_obssumm_filedb_200311.txt
    but only using the first two rows of data.
    """
    return {'filename': ['hsi_obssumm_20031101_139.fit',
                         'hsi_obssumm_20031102_144.fit'],
            'orb_st': [0, 0],
            'orb_end': [0, 0],
            'start_time': [datetime(2003, 11, 1, 0, 0),
                           datetime(2003, 11, 2, 0, 0)],
            'end_time': [datetime(2003, 11, 2, 0, 0),
                         datetime(2003, 11, 3, 0, 0)],
            'status_flag': [0, 0],
            'npackets': [0, 0]}


@mock.patch('sunpy.instr.rhessi.get_base_url',
            return_value='http://www.example.com')
@mock.patch('sunpy.instr.rhessi.parse_observing_summary_dbase_file',
            return_value=parsed_dbase())
@mock.patch('sunpy.instr.rhessi.get_observing_summary_dbase_file',
            return_value=('', {}))
def test_get_obssum_filename_one_day(mock_get_obssumm_dbase_file,
                                     mock_parse_obssumm_dbase_file,
                                     mock_get_base_url):
    """
    Given a time range of one day, make sure we get one days data back, i.e.
    one file.
    """
    tr = ('2003-11-01 01:00', '2003-11-01 02:00')
    filename = rhessi.get_observing_summary_filename(tr)

    assert len(filename) == 1
    assert filename[0] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031101_139.fits'


@mock.patch('sunpy.instr.rhessi.get_base_url',
            return_value='http://www.example.com')
@mock.patch('sunpy.instr.rhessi.parse_observing_summary_dbase_file',
            return_value=parsed_dbase())
@mock.patch('sunpy.instr.rhessi.get_observing_summary_dbase_file',
            return_value=('', {}))
def test_get_observing_summary_filename_two_days(mock_get_obssumm_dbase_file,
                                                 mock_parse_obssumm_dbase_file,
                                                 mock_get_base_url):
    """
    Given a time range of two days, make sure we get two files back, one
    for each day.
    """
    tr = ('2003-11-01', '2003-11-03')
    filenames = rhessi.get_observing_summary_filename(tr)

    assert len(filenames) == 2
    assert filenames[0] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031101_139.fits'
    assert filenames[1] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031102_144.fits'

# Test `rhessi.get_observing_summary_file(...)`


# Test `rhessi._build_energy_bands(...)`

@pytest.fixture
def raw_bands():
    """The RHESSI summary data standard energy bands."""
    return ['3 - 6', '6 - 12', '12 - 25', '25 - 50', '50 - 100', '100 - 300',
            '300 - 800', '800 - 7000', '7000 - 20000']


def test_build_energy_bands_no_match(raw_bands):
    """
    If an energy unit cannot be found in the `label` then raise
    a `ValueError`
    """
    with pytest.raises(ValueError):
        rhessi._build_energy_bands(label='Energy bands GHz', bands=raw_bands)


def test_build_energy_bands(raw_bands):
    """
    Success case.
    """
    built_ranges = rhessi._build_energy_bands(label='Energy bands (keV)',
                                              bands=raw_bands)

    assert built_ranges == ['3 - 6 keV', '6 - 12 keV', '12 - 25 keV',
                            '25 - 50 keV', '50 - 100 keV', '100 - 300 keV',
                            '300 - 800 keV', '800 - 7000 keV',
                            '7000 - 20000 keV']
