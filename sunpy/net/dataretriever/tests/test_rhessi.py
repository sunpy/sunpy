import socket
from datetime import datetime
from urllib.error import URLError
from urllib.request import urlretrieve

from unittest import mock
import pytest

import sunpy.instr.rhessi
import sunpy.net.dataretriever.sources.rhessi as rhessi
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import TimeRange, parse_time
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.dataretriever.client import QueryResponse

LCClient = rhessi.RHESSIClient()


@pytest.fixture
def one_day_timerange():
    return TimeRange(("2016/01/15", "2016/01/16"))


@pytest.fixture
def two_days_timerange():
    return TimeRange(("2016/01/15", "2016/01/17"))


def test_get_observing_summary_dbase_file_with_unsupported_start_time():
    """
    RHESSI summary files are not available for before 2002-02-01, ensure
    `ValueError` is raised.
    """
    with pytest.raises(ValueError):
        LCClient.get_observing_summary_dbase_file("2002/01/21")


@mock.patch('sunpy.net.dataretriever.sources.rhessi.urlretrieve', return_value=None)
@mock.patch('sunpy.net.dataretriever.sources.rhessi.get_base_url', return_value='http://www.example.com')
def test_get_observing_summary_dbase_file_build_correct_url(mock_get_base_url, mock_urlretrieve,
                                                            one_day_timerange):
    """
    This test ensures that we build the correct url which is then used
    to get the database file.
    """
    LCClient.get_observing_summary_dbase_file(one_day_timerange.start)
    mock_urlretrieve.assert_called_with(
        'http://www.example.com/dbase/hsi_obssumm_filedb_201601.txt')


# Test `rhessi.get_base_url()`

@mock.patch('sunpy.net.dataretriever.sources.rhessi.urlopen', return_value=None)
def test_get_base_url(mock_urlopen):
    """
    Success case, can successfully 'ping' first data_server
    """
    assert rhessi.get_base_url() == rhessi.data_servers[0]


@mock.patch('sunpy.net.dataretriever.sources.rhessi.urlopen', side_effect=URLError(''))
def test_get_base_url_on_urlerror(mock_urlopen):
    """
    If all tested URLs raise `URLError`, then raise an `IOError`
    """
    with pytest.raises(IOError):
        rhessi.get_base_url()


@mock.patch('sunpy.net.dataretriever.sources.rhessi.urlopen', side_effect=socket.timeout)
def test_get_base_url_on_timeout(mock_urlopen):
    """
    If all tested data servers timeout, then raise an `IOError`
    """
    with pytest.raises(IOError):
        rhessi.get_base_url()


# Test `rhessi.get_observing_summary_filename(...)`

def parsed_dbase():
    """
    The result of calling `parse_observing_summary_dbase_file(...)` on
    https://hesperia.gsfc.nasa.gov/hessidata/dbase/hsi_obssumm_filedb_200311.txt but
    only using the first two rows of data.
    """

    return {'filename': ['hsi_obssumm_20031101_139.fit',
                         'hsi_obssumm_20031102_144.fit',
                         'hsi_obssumm_20031103_148.fit',
                         'hsi_obssumm_20031104_150.fit',
                         'hsi_obssumm_20031105_138.fit'],
            'orb_st': [0, 0, 0, 0, 0],
            'orb_end': [0, 0, 0, 0, 0],
            'start_time': [datetime(2003, 11, 1, 0, 0),
                           datetime(2003, 11, 2, 0, 0),
                           datetime(2003, 11, 3, 0, 0),
                           datetime(2003, 11, 4, 0, 0),
                           datetime(2003, 11, 5, 0, 0)],
            'end_time': [datetime(2003, 11, 2, 0, 0),
                         datetime(2003, 11, 3, 0, 0),
                         datetime(2003, 11, 4, 0, 0),
                         datetime(2003, 11, 5, 0, 0),
                         datetime(2003, 11, 6, 0, 0)],
            'status_flag': [0, 0, 0, 0, 0],
            'npackets': [0, 0, 0, 0, 0]}


@pytest.mark.remote_data
def test_parsed_dbase():
    """
    Test that parsed_based still matches what is returned by rhessi.
    """
    filename, _ = urlretrieve("https://hesperia.gsfc.nasa.gov/hessidata/dbase/hsi_obssumm_filedb_200311.txt")
    dbase = sunpy.instr.rhessi.parse_observing_summary_dbase_file(filename)
    rows = {}
    for key in dbase.keys():
        rows[key] = dbase[key][:5]
    assert rows == parsed_dbase()

@mock.patch('sunpy.net.dataretriever.sources.rhessi.get_base_url', return_value='http://www.example.com')
@mock.patch('sunpy.instr.rhessi.parse_observing_summary_dbase_file', return_value=parsed_dbase())
@mock.patch('sunpy.net.dataretriever.sources.rhessi.RHESSIClient.get_observing_summary_dbase_file', return_value=('', {}))
def test_get_observing_summary_filename_one_day(mock_get_observing_summary_dbase_file,
                                                mock_parse_observing_summary_dbase_file,
                                                mock_get_base_url):
    """
    Given a time range of one day, make sure we get one days data back, i.e. one file.
    """
    filename = LCClient.get_observing_summary_filename(('2003-11-01', '2003-11-01T23:59:59'))

    assert len(filename) == 1
    assert filename[0] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031101_139.fits'


@mock.patch('sunpy.net.dataretriever.sources.rhessi.get_base_url', return_value='http://www.example.com')
@mock.patch('sunpy.instr.rhessi.parse_observing_summary_dbase_file', return_value=parsed_dbase())
@mock.patch('sunpy.net.dataretriever.sources.rhessi.RHESSIClient.get_observing_summary_dbase_file', return_value=('', {}))
def test_get_observing_summary_filename_two_days(mock_get_observing_summary_dbase_file,
                                                 mock_parse_observing_summary_dbase_file,
                                                 mock_get_base_url):
    """
    Given a time range of two days, make sure we get two files back, one
    for each day.
    """
    filenames = LCClient.get_observing_summary_filename(('2003-11-01', '2003-11-02T23:59:59'))

    assert len(filenames) == 2
    assert filenames[0] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031101_139.fits'
    assert filenames[1] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031102_144.fits'


def test_can_handle_query():
    ans1 = rhessi.RHESSIClient._can_handle_query(
        a.Time('2012/8/9', '2012/8/9'), a.Instrument('rhessi'))
    assert ans1 is True
    ans2 = rhessi.RHESSIClient._can_handle_query(a.Time('2013/2/7', '2013/2/7'))
    assert ans2 is False


@mock.patch('sunpy.net.dataretriever.sources.rhessi.get_base_url', return_value='http://www.example.com')
@mock.patch('sunpy.instr.rhessi.parse_observing_summary_dbase_file', return_value=parsed_dbase())
@mock.patch('sunpy.net.dataretriever.sources.rhessi.RHESSIClient.get_observing_summary_dbase_file', return_value=('', {}))
def test_query(mock_get_observing_summary_dbase_file,
               mock_parse_observing_summary_dbase_file,
               mock_get_base_url):
    qr1 = LCClient.search(a.Time('2003-11-01', '2003-11-03'), a.Instrument('rhessi'))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 3
    assert qr1.time_range().start.datetime == parse_time('2003/11/01')
    assert qr1.time_range().end.datetime == parse_time('2003/11/03T23:59:59.999')


@mock.patch('sunpy.net.dataretriever.sources.rhessi.get_base_url', return_value='http://www.example.com')
@mock.patch('sunpy.instr.rhessi.parse_observing_summary_dbase_file', return_value=parsed_dbase())
@mock.patch('sunpy.net.dataretriever.sources.rhessi.RHESSIClient.get_observing_summary_dbase_file', return_value=('', {}))
def test_fido_mock(mock_get_observing_summary_dbase_file,
                   mock_parse_observing_summary_dbase_file,
                   mock_get_base_url):
    qr = Fido.search(a.Time('2003-11-01', '2003-11-03'), a.Instrument('rhessi'))
    assert isinstance(qr, UnifiedResponse)
    assert qr._numfile == 3
