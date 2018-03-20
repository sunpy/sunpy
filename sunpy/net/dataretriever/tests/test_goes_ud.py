import datetime

import pytest
from hypothesis import given, example

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.tests.strategies import goes_time


@pytest.fixture
def LCClient():
    return goes.XRSClient()


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


@given(goes_time())
def test_can_handle_query(time):
    ans1 = goes.XRSClient._can_handle_query(time, Instrument('XRS'))
    assert ans1 is True
    ans2 = goes.XRSClient._can_handle_query(time)
    assert ans2 is False
    ans3 = goes.XRSClient._can_handle_query(time, Instrument('eve'))
    assert ans3 is False


def test_no_satellite(LCClient):
    with pytest.raises(ValueError):
        LCClient.search(Time("1950/01/01", "1950/02/02"), Instrument('XRS'))


def test_fixed_satellite(LCClient):
    ans1 = LCClient.search(a.Time("2017/01/01", "2017/01/02"),
                           a.Instrument('XRS'))

    for resp in ans1:
        assert "go15" in resp.url

    ans1 = LCClient.search(a.Time("2017/01/01", "2017/01/02"),
                           a.Instrument('XRS'),
                           a.goes.SatelliteNumber(13))

    for resp in ans1:
        assert "go13" in resp.url


@example(a.Time("2006-08-01", "2006-08-01"))
# This example tests a time range with a satellite jump and no overlap
@example(a.Time("2009-11-30", "2009-12-3"))
@given(goes_time())
def test_query(LCClient, time):
    qr1 = LCClient.search(time, Instrument('XRS'))
    assert isinstance(qr1, QueryResponse)
    # We only compare dates here as the start time of the qr will always be the
    # start of the day.
    assert qr1.time_range().start.date() == time.start.date()

    almost_day = datetime.timedelta(days=1, milliseconds=-1)
    end = datetime.datetime.combine(time.end.date(), datetime.time()) + almost_day
    assert qr1.time_range().end == end


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
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
def test_new_logic(LCClient):
    qr = LCClient.search(Time('2012/10/4', '2012/10/6'), Instrument('XRS'))
    res = LCClient.fetch(qr)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "time, instrument",
    [(a.Time("2012/10/4", "2012/10/6"), a.Instrument("goes")),
     (a.Time('2013/10/5', '2013/10/7'), a.Instrument("goes"))])
def test_fido(time, instrument):
    qr = Fido.search(a.Time('2012/10/4', '2012/10/6'), Instrument('XRS'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


@given(goes_time())
def test_time_for_url(LCClient, time):
    time = time.start.date().strftime("%Y/%m/%d")
    almost_day = datetime.timedelta(days=1, milliseconds=-1)

    tr = TimeRange(time, almost_day)
    url = LCClient._get_url_for_timerange(tr)
    times = LCClient._get_time_for_url(url)

    assert all([tr == t2 for t2 in times])
