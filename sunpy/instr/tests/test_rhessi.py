import platform
import textwrap
from unittest import mock
from distutils.version import LooseVersion

from os.path import sep
import numpy as np
import pytest

import sunpy.instr.rhessi as rhessi
import sunpy.io
import sunpy.map
from sunpy.data.test import get_test_filepath, rootdir
from sunpy.time import is_time_equal, parse_time


@pytest.fixture
def cross_month_timerange():
    """
    Time range which crosses a month boundary.

    Dbase files are monthly therefore this is to make sure that two
    dbase files are returned.
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
    Test that we get the observing summary database file with the content we
    expect.
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
    Ensure that all required data are extracted from the RHESSI observing
    summary database file mocked in ``hessi_data()``.
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
    """
    The RHESSI summary data standard energy bands.
    """
    return ['3 - 6', '6 - 12', '12 - 25', '25 - 50', '50 - 100', '100 - 300',
            '300 - 800', '800 - 7000', '7000 - 20000']


def test_build_energy_bands_no_match(raw_bands):
    """
    If an energy unit cannot be found in the ``label`` then raise a
    `ValueError`
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


def test_get_flare_list():
    src = rootdir if rootdir.endswith(sep) else (rootdir + sep)
    fa = 0.001  # float accuracy

    fl = rhessi.get_flare_list("2018-01-01", "2018-03-31", source=src)
    assert len(fl) == 24
    assert fl.loc[0]['ID_NUMBER'] == 18010601
    assert fl.loc[1]['START_TIME'].strftime("%Y-%m-%d %H:%M:%S") == "2018-01-18 07:29:24"
    assert fl.loc[2]['END_TIME'].strftime("%Y-%m-%d %H:%M:%S") == "2018-01-18 08:16:16"
    assert fl.loc[3]['PEAK_TIME'].strftime("%Y-%m-%d %H:%M:%S") == "2018-01-20 04:51:58"
    assert fl.loc[4]['BCK_TIME'][0].strftime("%Y-%m-%d %H:%M:%S") == "2018-01-20 04:52:28"
    assert fl.loc[5]['BCK_TIME'][1].strftime("%Y-%m-%d %H:%M:%S") == "2018-01-22 02:43:28"
    assert fl.loc[6]['IMAGE_TIME'][0].strftime("%Y-%m-%d %H:%M:%S") == "2018-02-04 21:14:12"
    assert fl.loc[7]['IMAGE_TIME'][1].strftime("%Y-%m-%d %H:%M:%S") == "2018-02-04 21:33:12"
    assert fl.loc[8]['ENERGY_RANGE_FOUND'][0] == pytest.approx(6.0, fa)
    assert fl.loc[9]['ENERGY_RANGE_FOUND'][1] == pytest.approx(12.0, fa)
    assert fl.loc[10]['ENERGY_HI'][0] == pytest.approx(6.0, fa)
    assert fl.loc[11]['ENERGY_HI'][1] == pytest.approx(50.0, fa)
    assert fl.loc[12]['PEAK_COUNTRATE'] == pytest.approx(6.0, fa)
    assert fl.loc[13]['BCK_COUNTRATE'] == pytest.approx(3.6386773586273193, fa)
    assert fl.loc[14]['TOTAL_COUNTS'] == pytest.approx(9808.0, fa)
    assert fl.loc[15]['PEAK_CORRECTION'] == pytest.approx(1.0, fa)
    assert fl.loc[16]['TOTAL_CORRECTION'] == pytest.approx(1.7999998331069946, fa)
    assert fl.loc[17]['POSITION'][0] == pytest.approx(-345.6, fa)
    assert fl.loc[18]['POSITION'][1] == 0
    assert fl.loc[19]['FILENAME'] == "hsi_20180209_152720_001.fits"
    assert fl.loc[20]['SFLAG1'] == 1
    assert fl.loc[21]['ACTIVE_REGION'] == 2699
    assert fl.loc[22]['GOES_CLASS'] == "B3.3"
    assert fl.loc[23]['ALT_ID'] == "20180303_040408"

    fl = rhessi.get_flare_list("2018-01-06 16:32:56", "2018-01-22 02:43:28", source=src)
    assert len(fl) == 6
    fl = rhessi.get_flare_list("2018-01-06 16:32:57", "2018-01-22 02:43:27", source=src)
    assert len(fl) == 4


def test_read_flare_list_file():
    src = (rootdir if rootdir.endswith(sep) else (rootdir + sep)) + "hessi_flare_list_201801.fits"
    fa = 0.001  # float accuracy
    fl = rhessi.read_flare_list_file(src)
    assert len(fl) == 6
    assert fl.loc[0]['PEAK_PHFLUX'] == pytest.approx(0.0, fa)
    assert fl.loc[1]['TOT_PHFLUENCE'] == pytest.approx(92538.421875, fa)
    assert fl.loc[2]['E_PHFLUENCE'] == pytest.approx(12945.400390625, fa)
    assert fl.loc[3]['PEAK_PHFLUX_SIGMA'] == pytest.approx(37513084.0, fa)
    assert fl.loc[4]['TOT_PHFLUENCE_SIGMA'] == pytest.approx(354753760.0, fa)
    assert fl.loc[5]['E_PHFLUENCE_SIGMA'] == pytest.approx(24626434.0, fa)


@pytest.mark.remote_data
def test_get_flare_list_remote():
    fl = rhessi.get_flare_list("2017-11-16", "2018-02-04")
    assert len(fl) == pytest.approx(13, 2)  # list may change in the future
