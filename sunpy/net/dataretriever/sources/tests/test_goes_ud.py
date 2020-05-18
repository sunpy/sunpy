import pytest
from hypothesis import given, settings

import astropy.units as u
from astropy.time import TimeDelta

import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net._attrs import Instrument, Time
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.tests.strategies import goes_time
from sunpy.time import is_time_equal, parse_time
from sunpy.time.timerange import TimeRange


@pytest.fixture
def LCClient():
    return goes.XRSClient()


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "timerange,url_start,url_end",
    [(TimeRange('1995/06/03', '1995/06/05'),
      'https://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits',
      'https://umbra.nascom.nasa.gov/goes/fits/1995/go07950605.fits'),
     (TimeRange('2008/06/02', '2008/06/04'),
      'https://umbra.nascom.nasa.gov/goes/fits/2008/go1020080602.fits',
      'https://umbra.nascom.nasa.gov/goes/fits/2008/go1020080604.fits')])
def test_get_url_for_time_range(LCClient, timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange, url_start, url_end",
                         [(TimeRange('1999/01/10', '1999/01/20'),
                           'https://umbra.nascom.nasa.gov/goes/fits/1999/go10990110.fits',
                           'https://umbra.nascom.nasa.gov/goes/fits/1999/go1019990120.fits')])
def test_get_overlap_urls(LCClient, timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert len(urls) == 9
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


@pytest.mark.remote_data
def test_no_satellite(LCClient):
    with pytest.raises(ValueError):
        LCClient.search(Time("1950/01/01", "1950/02/02"), Instrument('XRS'))


@pytest.mark.remote_data
def test_fixed_satellite(LCClient):
    ans1 = LCClient.search(a.Time("2017/01/01", "2017/01/02"),
                           a.Instrument.xrs)

    for resp in ans1:
        assert "go15" in resp.url

    ans1 = LCClient.search(a.Time("2017/01/01", "2017/01/02"),
                           a.Instrument.xrs,
                           a.goes.SatelliteNumber(13))

    for resp in ans1:
        assert "go13" in resp.url


@pytest.mark.parametrize("time", [
    Time('2005/4/27', '2005/4/27'),
    Time('2016/2/4', '2016/2/10')])
@pytest.mark.remote_data
def test_query(LCClient, time):
    qr1 = LCClient.search(time, Instrument('XRS'))
    assert isinstance(qr1, QueryResponse)
    # We only compare dates here as the start time of the qr will always be the
    # start of the day.
    assert qr1.time_range().start.strftime('%Y-%m-%d') == time.start.strftime('%Y-%m-%d')

    almost_day = TimeDelta(1*u.day - 1*u.millisecond)
    end = parse_time(time.end.strftime('%Y-%m-%d')) + almost_day
    assert is_time_equal(qr1.time_range().end, end)


@pytest.mark.remote_data
def test_query_error(LCClient):
    times = [a.Time("1983-05-01", "1983-05-02")]
    for time in times:
        with pytest.raises(ValueError):
            LCClient.search(time, Instrument('XRS'))


@pytest.mark.remote_data
@pytest.mark.parametrize("time, instrument", [
    (Time('1983/06/17', '1983/06/18'), Instrument('XRS')),
    (Time('2012/10/4', '2012/10/6'), Instrument('XRS')),
])
def test_get(LCClient, time, instrument):
    qr1 = LCClient.search(time, instrument)
    download_list = LCClient.fetch(qr1)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
def test_new_logic(LCClient):
    qr = LCClient.search(Time('2012/10/4', '2012/10/6'), Instrument('XRS'))
    download_list = LCClient.fetch(qr)
    assert len(download_list) == len(qr)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "time, instrument",
    [(a.Time("2012/10/4", "2012/10/5"), a.Instrument.goes)])
def test_fido(time, instrument):
    qr = Fido.search(time, Instrument('XRS'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


@settings(deadline=10000, max_examples=5)
@pytest.mark.remote_data
@given(goes_time())
def test_time_for_url(LCClient, time):
    time = time.start.strftime("%Y/%m/%d")
    almost_day = TimeDelta(1*u.day - 1*u.millisecond)

    tr = TimeRange(time, almost_day)
    url = LCClient._get_url_for_timerange(tr)
    times = LCClient._get_time_for_url(url)
    assert all([tr == t2 for t2 in times])


def test_attr_reg():
    assert a.Instrument.goes == a.Instrument("GOES")
    assert a.Instrument.xrs == a.Instrument("XRS")
    assert a.goes.SatelliteNumber.two == a.goes.SatelliteNumber("2")


def test_client_repr(LCClient):
    """
    Repr check
    """
    output = str(LCClient)
    assert output[:50] == 'sunpy.net.dataretriever.sources.goes.XRSClient\n\nPr'
