import pytest

from sunpy.time.timerange import TimeRange
from sunpy.time import parse_time
from sunpy.net._attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.gong_synoptic as gong_synoptic
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

from hypothesis import given
from sunpy.net.tests.strategies import range_time

GSClient = gong_synoptic.GongSynopticClient()


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange,url_start,url_end", [
    (TimeRange('2020/1/30', '2020/2/1'),
     'https://gong2.nso.edu/oQR/zqs/202001/mrzqs200130/mrzqs200130t0004c2227_349.fits.gz',
     'https://gong2.nso.edu/oQR/zqs/202001/mrzqs200131/mrzqs200131t2314c2227_323.fits.gz'
     ),
    (TimeRange('2020/4/21', '2020/4/22'),
     'https://gong2.nso.edu/oQR/zqs/202004/mrzqs200421/mrzqs200421t0004c2230_348.fits.gz',
     'https://gong2.nso.edu/oQR/zqs/202004/mrzqs200421/mrzqs200421t2314c2230_335.fits.gz'
     ),
    (TimeRange('2006/9/19', '2006/9/19 22:00'),
     'https://gong2.nso.edu/oQR/zqs/200609/mrzqs060919/mrzqs060919t1154c2048_323.fits.gz',
     'https://gong2.nso.edu/oQR/zqs/200609/mrzqs060919/mrzqs060919t1754c2048_320.fits.gz'
     )
])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = GSClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@given(range_time('2010-01-06'))
def test_can_handle_query(time):
    ans1 = gong_synoptic.GongSynopticClient._can_handle_query(
        time, Instrument('gong'))
    assert ans1 is True
    ans2 = gong_synoptic.GongSynopticClient._can_handle_query(time)
    assert ans2 is False
    ans3 = gong_synoptic.GongSynopticClient._can_handle_query(time, Instrument('goes'))
    assert ans3 is False


@pytest.mark.parametrize("time", [
    Time('2015/8/27', '2015/8/27 12:00'),
    Time('2020/2/4', '2020/2/5')])
@pytest.mark.remote_data
def test_query(time):
    qr1 = GSClient.search(time, Instrument('gong'))
    assert isinstance(qr1, QueryResponse)
    assert qr1.time_range().start == time.start
    assert qr1.time_range().end == time.end


@pytest.mark.remote_data
@pytest.mark.parametrize("time,instrument", [
    (Time('2013/8/27', '2013/8/27'), Instrument('gong')),
    (Time('2020/4/23 17:00', '2020/4/23 21:00'), Instrument('gong')),
])
def test_get(time, instrument):
    qr1 = GSClient.search(time, instrument)
    download_list = GSClient.fetch(qr1)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "time, instrument",
    [(a.Time('2019/10/4', '2019/10/4 2:00'), a.Instrument('gong')),
     (a.Time('2019/12/31 21:00', '2020/1/1'), a.Instrument('gong'))])
def test_fido(time, instrument):
    qr = Fido.search(time, instrument)
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
    for res in response:
        assert res.split('.')[-1] == 'fits'
