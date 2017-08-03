import pytest

from sunpy.time.timerange import TimeRange, parse_time
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

from hypothesis import given, example
from sunpy.net.tests.strategies import goes_time

LCClient = goes.GOESClient()


@pytest.mark.parametrize(
    "timerange,url_start,url_end",
    [(TimeRange('1995/06/03', '1995/06/05'),
      'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits',
      'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950605.fits'),
     (TimeRange('2008/06/02', '2008/06/04'),
      'http://umbra.nascom.nasa.gov/goes/fits/2008/go1020080602.fits',
      'http://umbra.nascom.nasa.gov/goes/fits/2008/go1020080604.fits')])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@given(goes_time())
def test_can_handle_query(time):
    ans1 = goes.GOESClient._can_handle_query(time, Instrument('goes'))
    assert ans1 is True
    ans2 = goes.GOESClient._can_handle_query(time)
    assert ans2 is False
    ans3 = goes.GOESClient._can_handle_query(time, Instrument('eve'))
    assert ans3 is False


def test_no_satellite():
    with pytest.raises(ValueError):
        LCClient.search(Time("1950/01/01", "1950/02/02"), Instrument('goes'))

@example(a.Time("2006-08-01", "2006-08-01"))
@example(a.Time("1983-05-01", "1983-05-02"))
# This example tests a time range with a satellite jump and no overlap
@example(a.Time("2009-11-30", "2009-12-3"))
@given(goes_time())
def test_query(time):
    tr = TimeRange(time.start, time.end)
    if parse_time("1983-05-01") in tr:
        with pytest.raises(ValueError):
            LCClient.search(time, Instrument('goes'))
    else:
        qr1 = LCClient.search(time, Instrument('goes'))
        assert isinstance(qr1, QueryResponse)
        assert qr1.time_range().start == time.start
        assert qr1.time_range().end == time.end


@pytest.mark.online
@pytest.mark.parametrize("time, instrument", [
    (Time('1983/06/17', '1983/06/18'), Instrument('goes')),
    (Time('2012/10/4', '2012/10/6'), Instrument('goes')),
])
def test_get(time, instrument):
    qr1 = LCClient.search(time, instrument)
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@pytest.mark.online
def test_new_logic():
    qr = LCClient.search(Time('2012/10/4', '2012/10/6'), Instrument('goes'))
    res = LCClient.fetch(qr)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr)


@pytest.mark.online
@pytest.mark.parametrize(
    "time, instrument",
    [(a.Time("2012/10/4", "2012/10/6"), a.Instrument("goes")),
     (a.Time('2013/10/5', '2013/10/7'), a.Instrument("goes"))])
def test_fido(time, instrument):
    qr = Fido.search(a.Time('2012/10/4', '2012/10/6'), Instrument('goes'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
