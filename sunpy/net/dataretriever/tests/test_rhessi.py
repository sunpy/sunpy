import pytest

import sunpy.net.dataretriever.sources.rhessi as rhessi
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse

LCClient = rhessi.RHESSIClient()


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange,url_start", [
    (TimeRange('2012/7/1', '2012/7/2'),
     'hessidata/metadata/catalog/hsi_obssumm_20120701_050.fits'),
    (TimeRange('2013/6/3', '2013/6/4'),
     'hessidata/metadata/catalog/hsi_obssumm_20130603_042.fits'),
    (TimeRange('2012/7/1', '2012/7/14'),
     'hessidata/metadata/catalog/hsi_obssumm_20120701_050.fits')
])
def test_get_url_for_time_range(timerange, url_start):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert url_start in urls[0]


def test_can_handle_query():
    ans1 = rhessi.RHESSIClient._can_handle_query(
        Time('2012/8/9', '2012/8/9'), Instrument('rhessi'))
    assert ans1 is True
    ans2 = rhessi.RHESSIClient._can_handle_query(Time('2013/2/7', '2013/2/7'))
    assert ans2 is False


@pytest.mark.remote_data
def test_query():
    qr1 = LCClient.search(Time('2011/4/9', '2011/4/10'), Instrument('rhessi'))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 1
    assert qr1.time_range().start == parse_time('2011/04/09')
    assert qr1.time_range().end == parse_time('2011/04/10')


@pytest.mark.remote_data
@pytest.mark.parametrize("time,instrument", [
    (Time('2012/11/27', '2012/11/28'), Instrument('rhessi')),
    (Time('2012/10/4', '2012/10/5'), Instrument('rhessi')),
])
def test_get(time, instrument):
    qr1 = LCClient.search(time, instrument)
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "time, instrument",
    [(a.Time('2012/10/4', '2012/10/6'), a.Instrument('rhessi')),
     (a.Time('2013/10/5', '2013/10/7'), a.Instrument('rhessi'))])
def test_fido(time, instrument):
    qr = Fido.search(time, instrument)
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
