import pytest
from hypothesis import given, example, settings
import tempfile

from astropy.time import TimeDelta
import astropy.units as u

from sunpy.time.timerange import TimeRange
from sunpy.net.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.tests.strategies import goes_time
from sunpy.time import parse_time, is_time_equal
from hypothesis import given
from sunpy.net.tests.strategies import time_attr


min_satellite_number = 16  # when SUVI was first included


sclient = goes.SUVIClient()


@given(time_attr())
def test_can_handle_query(time):
    ans1 = sclient._can_handle_query(time, a.Instrument('suvi'))
    assert ans1 is True
    ans2 = sclient._can_handle_query(time, a.Instrument('suvi'),
                                     a.Wavelength(131 * u.Angstrom))
    assert ans2 is True
    ans3 = sclient._can_handle_query(time, a.Instrument('suvi'),
                                     a.Wavelength(131 * u.Angstrom), a.Level(2))
    assert ans3 is True
    ans4 = sclient._can_handle_query(time)
    assert ans4 is False
    ans5 = sclient._can_handle_query(time, a.Instrument('aia'))
    assert ans5 is False
    ans6 = sclient._can_handle_query(time, a.Instrument('suvi'), SatelliteNumber(16))


def test_get_goes_sat_num():
    date = parse_time('2019/06/11 00:00')
    assert sclient._get_goes_sat_num(date) > min_satellite_number
    assert type(sclient._get_goes_sat_num(date)) is int


def test_get_goes_sat_num_error():
    date = parse_time('1800/06/11 00:00')
    with pytest.raises(ValueError):
        sclient._get_goes_sat_num(date)


def test_get_url_for_timerange_errors():
    """Check that unsupported values raise errors."""
    tr = TimeRange('2019/06/11 00:00', '2019/06/11 00:10')
    with pytest.raises(ValueError):
        sclient._get_url_for_timerange(tr, level=0)
    with pytest.raises(ValueError):
        sclient._get_url_for_timerange(tr, wavelength=100 * u.Angstrom)
    with pytest.raises(ValueError):
        sclient._get_url_for_timerange(tr, satellitenumber=1)


def mock_querry_object(start, end):
    """
    Creating a Query Response object and prefilling it with some information
    """
    # Creating a Query Response Object
    obj = {
        'TimeRange': TimeRange(parse_time(start), parse_time(end)),
        'Time_start': parse_time(start),
        'Time_end':  parse_time(end),
        'source': 'GOES',
        'instrument': 'SUVI',
        'physobs': 'flux',
        'provider': 'NOAA'
    }
    results = QueryResponse.create(obj, sclient._get_url_for_timerange(TimeRange(start, end)))
    results.client = sclient
    return results

start = '2019/05/25 00:50'
end = '2019/05/25 00:52'

@mock.patch('sunpy.net.dataretriever.sources.goes.SUVIClient.search',
            return_value=mock_querry_object(start, end))
def test_fido_query(mock_search):
    qr1 = sclient.search(a.Time(start, end), a.Instrument('suvi'),
                         a.Wavelength(171 * u.Angstrom))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 1
    assert qr1.time_range().start == parse_time(start)
    assert qr1.time_range().end == parse_time(end)


@pytest.mark.remote_data
def test_fetch_working():
    """
    Tests if the online server for fermi_gbm is working.
    This also checks if the mock is working well.
    """
    start = '2019/05/25 00:50'
    end = '2019/05/25 00:52'
    qr1 = sclient.search(a.Time(start, end), a.Instrument('suvi'))

    # Mock QueryResponse object
    mock_qr = mock_querry_object(start, end)

    # Compare if two objects have the same attribute

    mock_qr = mock_qr[0]
    qr = qr1[0]

    assert mock_qr.source == qr.source
    assert mock_qr.provider == qr.provider
    assert mock_qr.physobs == qr.physobs
    assert mock_qr.instrument == qr.instrument
    assert mock_qr.url == qr.url
    assert mock_qr.time == qr.time

    assert qr1.time_range() == TimeRange(start, end)

    download_list = sclient.fetch(qr1, path=tempfile.TemporaryDirectory())
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize("start, end,wave,expected_num_files",
                         ('2019/05/25 00:50', '2019/05/25 00:52', 94, 1),
                         ('2019/05/25 00:50', '2019/05/25 00:52', 131, 1),
                         ('2019/05/25 00:50', '2019/05/25 00:52', 171, 1),
                         ('2019/05/25 00:50', '2019/05/25 00:52', 195, 1),
                         ('2019/05/25 00:50', '2019/05/25 00:52', 284, 1),
                         ('2019/05/25 00:50', '2019/05/25 00:52', 304, 1))
def test_get_url_for_time_range_level2(start, end, wave, expected_num_files):
    urls = sclient._get_url_for_timerange(TimeRange(start, end),
                                          wavelength=wave * u.Angstrom,
                                          level=2)
    assert isinstance(urls, list)
    assert len(urls) == expected_num_files

@pytest.mark.remote_data
@pytest.mark.parametrize("start, end,wave,expected_num_files",
                         ('2019/05/25 00:50', '2019/05/25 00:52', 6)
                         )
def test_get_url_for_time_range_level2_allwave(start, end, expected_num_files):
    """check that we get all wavelengths if no wavelength is given"""
    urls = sclient._get_url_for_timerange(TimeRange(start, end), level=2)
    assert isinstance(urls, list)
    assert len(urls) == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize("start, end,wave,expected_num_files",
                         ('2019/05/25 00:50', '2019/05/25 00:52', 94, 6),
                         ('2019/05/25 00:50', '2019/05/25 00:52', 131, 6),
                         ('2019/05/25 00:50', '2019/05/25 00:52', 171, 2),
                         ('2019/05/25 00:50', '2019/05/25 00:52', 195, 9),
                         ('2019/05/25 00:50', '2019/05/25 00:52', 284, 2),
                         ('2019/05/25 00:50', '2019/05/25 00:52', 304, 5))
def test_get_url_for_time_range_level1b(start, end, wave, expected_num_files):
    """check that we get all wavelengths if no wavelength is given"""
    urls = sclient._get_url_for_timerange(TimeRange(start, end), level='1b')
    assert isinstance(urls, list)
    assert len(urls) == expected_num_files
