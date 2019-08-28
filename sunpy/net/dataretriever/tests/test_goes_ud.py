import pytest
from hypothesis import given, example, settings

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
from hypothesis import given
from sunpy.net.tests.strategies import time_attr


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


@settings(deadline=50000)
@example(a.Time("2006-08-01", "2006-08-01"))
# This example tests a time range with a satellite jump and no overlap
@example(a.Time("2009-11-30", "2009-12-3"))
@given(goes_time())
def test_query(LCClient, time):
    qr1 = LCClient.search(time, Instrument('XRS'))
    assert isinstance(qr1, QueryResponse)
    # We only compare dates here as the start time of the qr will always be the
    # start of the day.
    assert qr1.time_range().start.strftime('%Y-%m-%d') == time.start.strftime('%Y-%m-%d')

    almost_day = TimeDelta(1*u.day - 1*u.millisecond)
    end = parse_time(time.end.strftime('%Y-%m-%d')) + almost_day
    assert is_time_equal(qr1.time_range().end, end)


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
    [(a.Time("2012/10/4", "2012/10/6"), a.Instrument("goes")),
     (a.Time('2013/10/5', '2013/10/7'), a.Instrument("goes"))])
def test_fido(time, instrument):
    qr = Fido.search(a.Time('2012/10/4', '2012/10/6'), Instrument('XRS'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile

@settings(deadline=50000)
@given(goes_time())
def test_time_for_url(LCClient, time):
    time = time.start.strftime("%Y/%m/%d")
    almost_day = TimeDelta(1*u.day - 1*u.millisecond)

    tr = TimeRange(time, almost_day)
    url = LCClient._get_url_for_timerange(tr)
    times = LCClient._get_time_for_url(url)

    assert all([tr == t2 for t2 in times])


sclient = goes.SUVIClient()


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange,url_start,url_end",
                         [(TimeRange('2019/05/13 00:00', '2019/05/13 01:00'),
                          'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/suvi-l2-ci171/2019/05/13/dr_suvi-l2-ci171_g16_s20190513T000000Z_e20190513T000400Z_v1-0-0.fits',
                           'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/suvi-l2-ci171/2019/05/13/dr_suvi-l2-ci171_g16_s20190513T005600Z_e20190513T010000Z_v1-0-0.fits'),
                          (TimeRange('2019/05/13 00:00', '2019/05/13 00:00'),
                           'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/suvi-l2-ci171/2019/06/11/dr_suvi-l2-ci171_g16_s20190611T000000Z_e20190611T000400Z_v1-0-0.fits',
                           'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/suvi-l2-ci171/2019/06/11/dr_suvi-l2-ci171_g16_s20190611T000800Z_e20190611T001200Z_v1-0-0.fits')])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = sclient._get_url_for_timerange(timerange, wavelength=171 * u.Angstrom, level=2)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@given(time_attr())
def test_can_handle_query(time):
    ans1 = sclient._can_handle_query(time, a.Instrument('suvi'))
    assert ans1 is True
    ans2 = sclient._can_handle_query(time, a.Instrument('suvi'),
                                      a.Wavelength(131 * u.Angstrom))
    assert ans2 is True

    ans4 = sclient._can_handle_query(time)
    assert ans4 is False


@pytest.mark.remote_data
@pytest.mark.parametrize("time,instrument", [
    (a.Time('2019/05/13 00:00', '2019/05/13 00:10'), a.Instrument('suvi')),
])
def test_query(time, instrument):
    qr1 = sclient.search(time, instrument)
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 3
    assert qr1.time_range().start == time.start
    assert qr1.time_range().end == time.end


@pytest.mark.remote_data
@pytest.mark.parametrize("time,instrument", [
    (a.Time('2019/05/13 00:00', '2019/05/13 01:00'), a.Instrument('suvi')),
])
def test_get(time, instrument):
    qr1 = sclient.search(time, instrument)
    download_list = sclient.fetch(qr1)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    'query',
    [(a.Time('2019/05/13 00:00', '2019/05/13 01:00') & a.Instrument('suvi') & a.Wavelength(195 * u.Angstrom))])
def test_fido(query):
    qr = Fido.search(query)
    client = qr.get_response(0).client
    assert isinstance(qr, UnifiedResponse)
    assert type(client) == type(sclient)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
