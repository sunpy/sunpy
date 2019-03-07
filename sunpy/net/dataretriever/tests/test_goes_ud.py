import pytest
from unittest import mock
from datetime import datetime, timedelta

from hypothesis import given, example, settings

from astropy.tests.helper import assert_quantity_allclose
from astropy.time import TimeDelta
import astropy.units as u

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.tests.strategies import goes_time
from sunpy.time import parse_time, is_time_equal


LCClient = goes.XRSClient()


def mock_query_object(start_date, end_date):
    """
    Creation of a QueryResponse object, and prefill some
    downloaded data from goes.XRSClient().fetch(Time('20 ..)
    """
    # Create a mock QueryResponse object
    map_ = {
        'TimeRange' : TimeRange(parse_time(start_date), parse_time(end_date)),
        'Time_start': parse_time(start_date),
        'Time_end':  parse_time(end_date),
        'source': 'nasa',
        'instrument': 'goes',
        'physobs': 'irradiance',
        'provider': 'sdac'
    }

    st_datetime = datetime.strptime(start_date, '%Y/%m/%d')
    ed_datetime = datetime.strptime(end_date, '%Y/%m/%d') + timedelta(days=1)
    time_range = TimeRange(st_datetime, ed_datetime).split(int((ed_datetime-st_datetime).days))

    resp = QueryResponse.create(map_,
                                LCClient._get_url_for_timerange(
                                    map_['TimeRange']),
                                time=time_range)
    # Attach the client with the QueryResponse
    resp.client = LCClient
    return resp


@pytest.mark.remote_data
def test_fetch_working():
    """
    Tests if the online server for noaa is working.
    Uses the url : ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt
    """
    qr1 = LCClient.search(Time('1995/06/03', '1995/06/05'),
                          Instrument('XRS'))

    # Mock QueryResponse object
    mock_qr = mock_query_object('1995/06/03', '1995/06/05')

    # Compare if two objects have the same attribute

    mock_qr = mock_qr[0]
    qr = qr1[0]

    assert mock_qr.source == qr.source
    assert mock_qr.provider == qr.provider
    assert mock_qr.physobs == qr.physobs
    assert mock_qr.instrument == qr.instrument
    assert mock_qr.url == qr.url

    # Here the mock object contains `datetime.datetime` object by default,
    # so no parsing required

    sdate1 = mock_qr.time.start.value
    sdate2 = datetime.strptime(qr.time.start.value, "%Y-%m-%dT%H:%M:%S.%f")

    edate1 = mock_qr.time.end.value
    edate2 = datetime.strptime(qr.time.end.value, "%Y-%m-%dT%H:%M:%S.%f")

    assert abs(sdate2 - sdate1) < timedelta(seconds=1)
    assert abs(edate2 - edate1) < timedelta(seconds=1)


@pytest.mark.parametrize(
    "timerange,url_start,url_end",
    [(TimeRange('1995/06/03', '1995/06/05'),
      'https://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits',
      'https://umbra.nascom.nasa.gov/goes/fits/1995/go07950605.fits'),
     (TimeRange('2008/06/02', '2008/06/04'),
      'https://umbra.nascom.nasa.gov/goes/fits/2008/go1020080602.fits',
      'https://umbra.nascom.nasa.gov/goes/fits/2008/go1020080604.fits')])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@given(goes_time())
def test_can_handle_query(time):
    ans1 = goes.XRSClient._can_handle_query(time, Instrument('XRS'))
    assert ans1 is True
    ans2 = goes.XRSClient._can_handle_query(time)
    assert ans2 is False
    ans3 = goes.XRSClient._can_handle_query(time, Instrument('eve'))
    assert ans3 is False


def test_no_satellite():
    with pytest.raises(ValueError):
        LCClient.search(Time("1950/01/01", "1950/02/02"), Instrument('XRS'))


def test_fixed_satellite():
    ans1 = LCClient.search(a.Time("2017/01/01", "2017/01/02"),
                           a.Instrument('XRS'))

    for resp in ans1:
        assert "go15" in resp.url

    ans1 = LCClient.search(a.Time("2017/01/01", "2017/01/02"),
                           a.Instrument('XRS'),
                           a.goes.SatelliteNumber(13))

    for resp in ans1:
        assert "go13" in resp.url


@settings(deadline=50000)
@example(a.Time("2006-08-01", "2006-08-01"))
# This example tests a time range with a satellite jump and no overlap
@example(a.Time("2009-11-30", "2009-12-3"))
@given(goes_time())
def test_query(time):
    qr1 = LCClient.search(time, Instrument('XRS'))
    assert isinstance(qr1, QueryResponse)
    # We only compare dates here as the start time of the qr will always be the
    # start of the day.
    assert qr1.time_range().start.strftime('%Y-%m-%d') == time.start.strftime('%Y-%m-%d')

    almost_day = TimeDelta(1*u.day - 1*u.millisecond)
    end = parse_time(time.end.strftime('%Y-%m-%d')) + almost_day
    assert is_time_equal(qr1.time_range().end, end)


def test_query_error():
    times = [a.Time("1983-05-01", "1983-05-02")]
    for time in times:
        with pytest.raises(ValueError):
            LCClient.search(time, Instrument('XRS'))


@mock.patch('sunpy.net.fido_factory.Fido.fetch',
            side_effect=(UnifiedResponse(
                mock_query_object('1983/06/17', '1983/06/18'))))
@mock.patch('sunpy.net.download.Results.wait', return_value=['/home/yash/sunpy/data/go06830617.fits',
                                                     '/home/yash/sunpy/data/go06830618.fits'])
def test_get(mock_result, mock_fetch):
    qr1 = LCClient.search(Time('1983/06/17', '1983/06/18'), Instrument('XRS'))
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@mock.patch('sunpy.net.fido_factory.Fido.fetch',
            side_effect=(UnifiedResponse(
                mock_query_object('2012/10/4', '2012/10/6'))))
@mock.patch('sunpy.net.download.Results.wait', return_value=['/home/yash/sunpy/data/go1520121004.fits',
 '/home/yash/sunpy/data/go1520121006.fits', '/home/yash/sunpy/data/go1520121005.fits'])
def test_get(mock_result, mock_fetch):
    qr1 = LCClient.search(Time('2012/10/4', '2012/10/6'), Instrument('XRS'))
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@mock.patch('sunpy.net.fido_factory.Fido.fetch',
            side_effect=(UnifiedResponse(
                mock_query_object('2012/10/4', '2012/10/6'))))
@mock.patch('sunpy.net.download.Results.wait', return_value=['/home/yash/sunpy/data/go1520121004.fits',
 '/home/yash/sunpy/data/go1520121006.fits', '/home/yash/sunpy/data/go1520121005.fits'])
def test_new_logic(mock_result, mock_fetch):
    qr = LCClient.search(Time('2012/10/4', '2012/10/6'), Instrument('XRS'))
    res = LCClient.fetch(qr)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr)



@mock.patch('sunpy.net.fido_factory.Fido.fetch',
            side_effect=(UnifiedResponse(
                mock_query_object('2012/10/4', '2012/10/6'))))
@mock.patch('sunpy.net.download.Results.wait', return_value=['/home/yash/sunpy/data/go1520121006.fits',
'/home/yash/sunpy/data/go1520121004.fits', '/home/yash/sunpy/data/go1520121005.fits']
)
def test_fido(time, instrument):
    qr = Fido.search(a.Time('2012/10/4', '2012/10/6'), Instrument('XRS'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


@mock.patch('sunpy.net.fido_factory.Fido.fetch',
            side_effect=(UnifiedResponse(
                mock_query_object('2013/10/5', '2013/10/7'))))
@mock.patch('sunpy.net.download.Results.wait', return_value=['/home/yash/sunpy/data/go1520131007.fits',
'/home/yash/sunpy/data/go1520131005.fits', '/home/yash/sunpy/data/go1520131006.fits']
)
def test_fido(time, instrument):
    qr = Fido.search(a.Time('2013/10/5', '2013/10/7'), Instrument('XRS'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


@settings(deadline=50000)
@given(goes_time())
def test_time_for_url(time):
    time = time.start.strftime("%Y/%m/%d")
    almost_day = TimeDelta(1*u.day - 1*u.millisecond)

    tr = TimeRange(time, almost_day)
    url = LCClient._get_url_for_timerange(tr)
    times = LCClient._get_time_for_url(url)

    assert all([tr == t2 for t2 in times])
