# -*- coding: utf-8 -*-
"""
Unit tests for `sunpy.instr.rhessi`
"""
import platform
import textwrap
from distutils.version import LooseVersion

from unittest import mock
import numpy as np
import pytest

import sunpy.io
import sunpy.map
from sunpy.data.test import get_test_filepath
import sunpy.instr.rhessi as rhessi
from sunpy.time import parse_time, is_time_equal


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
    amap = rhessi.backprojection(get_test_filepath(test_filename))
    assert isinstance(amap, sunpy.map.GenericMap)
    assert is_time_equal(amap.date, parse_time((2002, 2, 20, 11, 6, 21)))


def test_parse_obssum_dbase_file():
    fname = get_test_filepath("hsi_obssumm_filedb_201104.txt")
    obssum = rhessi.parse_observing_summary_dbase_file(fname)
    assert obssum['filename'][0] == 'hsi_obssumm_20110401_043.fit'
    assert obssum['filename'][-1] == 'hsi_obssumm_20110430_029.fit'

    assert obssum['orb_st'][0] == 0
    assert obssum['orb_st'][-1] == 0

    assert obssum['orb_end'][0] == 0
    assert obssum['orb_end'][-1] == 0

    assert obssum['start_time'][0] == parse_time((2011, 4, 1, 0, 0, 0))
    assert obssum['start_time'][-1] == parse_time((2011, 4, 30, 0, 0, 0))

    assert obssum['end_time'][0] == parse_time((2011, 4, 2, 0, 0, 0))
    assert obssum['end_time'][-1] == parse_time((2011, 5, 1, 0, 0, 0))

    assert obssum['status_flag'][0] == 0
    assert obssum['status_flag'][-1] == 0

    assert obssum['npackets'][0] == 0
    assert obssum['npackets'][-1] == 0


def test_parse_observing_summary_dbase_file():
    """
    Test that we get the observing summary dbase file with the content
    we expect.
    """
    obssum = rhessi.parse_observing_summary_dbase_file(get_test_filepath("hsi_obssumm_filedb_201104.txt"))

    assert obssum['filename'][0][0:20] == 'hsi_obssumm_20110401'
    assert obssum['filename'][1][0:20] == 'hsi_obssumm_20110402'

    assert obssum['orb_st'][0] == 0
    assert obssum['orb_st'][-1] == 0

    assert obssum['orb_end'][0] == 0
    assert obssum['orb_end'][-1] == 0

    assert obssum['start_time'][0] == parse_time((2011, 4, 1, 0, 0, 0))
    assert obssum['start_time'][-1] == parse_time((2011, 4, 30, 0, 0, 0))

    assert obssum['end_time'][0] == parse_time((2011, 4, 2, 0, 0, 0))
    assert obssum['end_time'][-1] == parse_time((2011, 5, 1, 0, 0, 0))

    assert obssum['status_flag'][0] == 0
    assert obssum['status_flag'][-1] == 0

    assert obssum['npackets'][0] == 0
    assert obssum['npackets'][-1] == 0


def test_get_parse_obssum_hdulist():
    hdulist = sunpy.io.read_file(get_test_filepath('hsi_obssumm_20110404_042.fits.gz'))
    header, _data = rhessi.parse_observing_summary_hdulist(hdulist)
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


# Test `rhessi.parse_obssumm_dbase_file(...)`


def hessi_data():
    return textwrap.dedent("""\
         HESSI Filedb File:
         Created: 1972-04-14T12:41:26.000
         Number of Files:           2
                             Filename  Orb_st Orb_end         Start_time           End_time Status_flag    Npackets Drift_start   Drift_end Data source
         hsi_obssumm_19721101_139.fit       7       8 01-Nov-72 00:00:00 02-Nov-72 00:00:00           3           2       0.000       0.000
         hsi_obssumm_19721102_144.fit       9      10 02-Nov-72 00:00:00 03-Nov-72 00:00:00           4           1       0.000       0.000
         """)


def test_parse_observing_summary_dbase_file_mock():
    """
    Ensure that all required data are extracted from the RHESSI
    observing summary database file mocked in `hessi_data()`
    """
    # We need to mock this test differently for <= 3.7.0 and below.
    if LooseVersion(platform.python_version()) <= LooseVersion("3.7.0"):
        mock_file = mock.mock_open()
        mock_file.return_value.__iter__.return_value = hessi_data().splitlines()
    else:
        mock_file = mock.mock_open(read_data=hessi_data())

    dbase_data = {}
    with mock.patch('sunpy.instr.rhessi.open', mock_file, create=True):
        dbase_data = rhessi.parse_observing_summary_dbase_file(None)

    assert len(dbase_data.keys()) == 7

    # verify each of the 7 fields
    assert dbase_data['filename'] == ['hsi_obssumm_19721101_139.fit',
                                      'hsi_obssumm_19721102_144.fit']
    assert dbase_data['orb_st'] == [7, 9]
    assert dbase_data['orb_end'] == [8, 10]
    assert dbase_data['start_time'] == [parse_time((1972, 11, 1, 0, 0)),
                                        parse_time((1972, 11, 2, 0, 0))]
    assert dbase_data['end_time'] == [parse_time((1972, 11, 2, 0, 0)),
                                      parse_time((1972, 11, 3, 0, 0))]
    assert dbase_data['status_flag'] == [3, 4]
    assert dbase_data['npackets'] == [2, 1]


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
