import re
import socket
from unittest import mock
from http.client import RemoteDisconnected
from urllib.error import URLError

import pytest

import sunpy.net.dataretriever.sources.rhessi as rhessi
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.tests.helpers import no_vso
from sunpy.time import TimeRange, parse_time


@pytest.fixture
def LCClient():
    return rhessi.RHESSIClient()


@pytest.fixture
def one_day_timerange():
    return TimeRange(("2016/01/15", "2016/01/16"))


@pytest.fixture
def two_days_timerange():
    return TimeRange(("2016/01/15", "2016/01/17"))


def test_get_observing_summary_dbase_file_with_unsupported_start_time(LCClient):
    """
    RHESSI summary files are not available for before 2002-02-01, ensure
    `ValueError` is raised.
    """
    with pytest.raises(ValueError, match="RHESSI summary files are not available before 2002-02-01"):
        LCClient.get_observing_summary_dbase_file("2002/01/21")


@mock.patch('sunpy.net.dataretriever.sources.rhessi.urlretrieve', return_value=None)
@mock.patch('sunpy.net.dataretriever.sources.rhessi.get_base_url', return_value='http://www.example.com/')
def test_get_observing_summary_dbase_file_build_correct_url(mock_get_base_url, mock_urlretrieve,
                                                            one_day_timerange, LCClient):
    """
    This test ensures that we build the correct url which is then used
    to get the database file.
    """
    LCClient.get_observing_summary_dbase_file(one_day_timerange.start)
    mock_urlretrieve.assert_called_with(
        'http://www.example.com/dbase/hsi_obssumm_filedb_201601.txt')


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
    with pytest.raises(OSError, match=re.escape("Unable to find an online HESSI server from")):
        rhessi.get_base_url()


@mock.patch('sunpy.net.dataretriever.sources.rhessi.urlopen', side_effect=socket.timeout)
def test_get_base_url_on_timeout(mock_urlopen):
    """
    If all tested data servers timeout, then raise an `IOError`
    """
    with pytest.raises(OSError, match=re.escape("Unable to find an online HESSI server from ('https://hesperia.gsfc.nasa.gov/hessidata/', 'http://hessi.ssl.berkeley.edu/hessidata/', 'http://soleil.i4ds.ch/hessidata/')")):
        rhessi.get_base_url()


@mock.patch('sunpy.net.dataretriever.sources.rhessi.urlopen', side_effect=RemoteDisconnected(''))
def test_get_base_url_on_remote_disconnected(mock_urlopen):
    """
    If all tested data servers timeout, then raise an `IOError`
    """
    with pytest.raises(OSError, match=re.escape("Unable to find an online HESSI server from ('https://hesperia.gsfc.nasa.gov/hessidata/', 'http://hessi.ssl.berkeley.edu/hessidata/', 'http://soleil.i4ds.ch/hessidata/')")):
        rhessi.get_base_url()


def parsed_dbase():
    """
    The result of calling `parse_observing_summary_dbase_file(...)` on
    https://hesperia.gsfc.nasa.gov/hessidata/dbase/hsi_obssumm_filedb_200311.txt but
    only using the first two rows of data.
    """

    return {'filename': ['hsi_obssumm_20031101_157.fit',
                         'hsi_obssumm_20031102_161.fit',
                         'hsi_obssumm_20031103_165.fit',
                         'hsi_obssumm_20031104_168.fit',
                         'hsi_obssumm_20031105_156.fit'],
            'orb_st': [65536, 65536, 65536, 65536, 65536],
            'orb_end': [65536, 65536, 65536, 65536, 65536],
            'start_time': [parse_time('2003-11-01T00:00:00.000'),
                           parse_time('2003-11-02T00:00:00.000'),
                           parse_time('2003-11-03T00:00:00.000'),
                           parse_time('2003-11-04T00:00:00.000'),
                           parse_time('2003-11-05T00:00:00.000')],
            'end_time': [parse_time('2003-11-02T00:00:00.000'),
                         parse_time('2003-11-03T00:00:00.000'),
                         parse_time('2003-11-04T00:00:00.000'),
                         parse_time('2003-11-05T00:00:00.000'),
                         parse_time('2003-11-06T00:00:00.000')],
            'status_flag': [0, 0, 0, 0, 0],
            'npackets': [0, 0, 0, 0, 0]}


@mock.patch('sunpy.net.dataretriever.sources.rhessi.get_base_url', return_value='http://www.example.com/')
@mock.patch('sunpy.net.dataretriever.sources.rhessi.parse_observing_summary_dbase_file', return_value=parsed_dbase())
@mock.patch('sunpy.net.dataretriever.sources.rhessi.RHESSIClient.get_observing_summary_dbase_file', return_value=('', {}))
def test_get_observing_summary_filename_one_day(mock_get_observing_summary_dbase_file,
                                                mock_parse_observing_summary_dbase_file,
                                                mock_get_base_url, LCClient):
    """
    Given a time range of one day, make sure we get one days data back, i.e. one file.
    """
    filename = LCClient.get_observing_summary_filename(('2003-11-01', '2003-11-01T23:59:59'))
    assert len(filename) == 1
    assert filename[0] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031101_157.fits'


@mock.patch('sunpy.net.dataretriever.sources.rhessi.get_base_url', return_value='http://www.example.com/')
@mock.patch('sunpy.net.dataretriever.sources.rhessi.parse_observing_summary_dbase_file', return_value=parsed_dbase())
@mock.patch('sunpy.net.dataretriever.sources.rhessi.RHESSIClient.get_observing_summary_dbase_file', return_value=('', {}))
def test_get_observing_summary_filename_two_days(mock_get_observing_summary_dbase_file,
                                                 mock_parse_observing_summary_dbase_file,
                                                 mock_get_base_url, LCClient):
    """
    Given a time range of two days, make sure we get two files back, one
    for each day.
    """
    filenames = LCClient.get_observing_summary_filename(('2003-11-01', '2003-11-02T23:59:59'))

    assert len(filenames) == 2
    assert filenames[0] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031101_157.fits'
    assert filenames[1] == 'http://www.example.com/metadata/catalog/hsi_obssumm_20031102_161.fits'


def test_can_handle_query(LCClient):
    ans1 = LCClient._can_handle_query(
        a.Time('2012/8/9', '2012/8/9'), a.Instrument.rhessi)
    assert ans1 is True
    ans2 = LCClient._can_handle_query(a.Time('2013/2/7', '2013/2/7'))
    assert ans2 is False


@mock.patch('sunpy.net.dataretriever.sources.rhessi.get_base_url', return_value='http://www.example.com/')
@mock.patch('sunpy.net.dataretriever.sources.rhessi.parse_observing_summary_dbase_file', return_value=parsed_dbase())
@mock.patch('sunpy.net.dataretriever.sources.rhessi.RHESSIClient.get_observing_summary_dbase_file', return_value=('', {}))
def test_query(mock_get_observing_summary_dbase_file,
               mock_parse_observing_summary_dbase_file,
               mock_get_base_url, LCClient):
    qr1 = LCClient.search(a.Time('2003-11-01', '2003-11-03'), a.Instrument.rhessi)
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 3
    assert qr1.time_range().start.datetime == parse_time('2003/11/01').datetime
    assert qr1.time_range().end.datetime == parse_time('2003/11/03T23:59:59.999').datetime


@no_vso
@mock.patch('sunpy.net.dataretriever.sources.rhessi.get_base_url', return_value='http://www.example.com/')
@mock.patch('sunpy.net.dataretriever.sources.rhessi.parse_observing_summary_dbase_file', return_value=parsed_dbase())
@mock.patch('sunpy.net.dataretriever.sources.rhessi.RHESSIClient.get_observing_summary_dbase_file', return_value=('', {}))
def test_fido_mock(mock_get_observing_summary_dbase_file,
                   mock_parse_observing_summary_dbase_file,
                   mock_get_base_url):
    qr = Fido.search(a.Time('2003-11-01', '2003-11-03'), a.Instrument.rhessi)
    assert isinstance(qr, UnifiedResponse)
    assert qr._numfile == 3


def test_attr_reg():
    assert a.Instrument.rhessi == a.Instrument('RHESSI')
    assert a.Physobs.summary_lightcurve == a.Physobs("summary_lightcurve")


def test_client_repr(LCClient):
    """
    Repr check
    """
    output = str(LCClient)
    assert output[:51] == 'sunpy.net.dataretriever.sources.rhessi.RHESSIClient'


def mock_query_object(LCClient):
    """
    Creating a Query Response object and prefilling it with some information
    """
    start = '2016/1/1'
    end = '2016/1/1 23:59:59'
    obj = {
        'Start Time': parse_time(start),
        'End Time': parse_time(end),
        'Instrument': 'RHESSI',
        'Physobs': 'irradiance',
        'Source': 'RHESSI',
        'Provider': 'NASA',
        'url': ('https://hesperia.gsfc.nasa.gov/hessidata/metadata/'
                'catalog/hsi_obssumm_20160101_078.fits')
    }
    results = QueryResponse([obj], client=LCClient)
    return results


def test_show(LCClient):
    mock_qr = mock_query_object(LCClient)
    qrshow0 = mock_qr.show()
    qrshow1 = mock_qr.show('Start Time', 'Instrument')
    allcols = {'Start Time', 'End Time', 'Instrument', 'Physobs', 'Source', 'Provider', 'url'}
    assert not allcols.difference(qrshow0.colnames)
    assert qrshow1.colnames == ['Start Time', 'Instrument']
    assert qrshow0['Instrument'][0] == 'RHESSI'
