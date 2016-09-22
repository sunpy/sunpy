import datetime
import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument, Source
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.norh as norh
from sunpy.net.dataretriever.downloader_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

from hypothesis import given
from .strategies import time_attr

LCClient = norh.NoRHClient()


@pytest.mark.parametrize("timerange,url_start,url_end", [
    (TimeRange('2012/4/21', '2012/4/21'),
     'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/04/tca120421',
     'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/04/tca120421'
     ),
    (TimeRange('2012/12/1', '2012/12/2'),
     'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/12/tca121201',
     'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/12/tca121202'
     ),
    (TimeRange('2012/3/7', '2012/3/14'),
     'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/03/tca120307',
     'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/03/tca120314'
     )
])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


def test_get_url_for_date():
    url = LCClient._get_url_for_date(datetime.date(2011, 3, 14))
    assert url == 'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2011/03/tca110314'


@given(time_attr())
def test_can_handle_query(time):
    ans1 = norh.NoRHClient._can_handle_query(time, Instrument('norh'))
    assert ans1 is True
    ans2 = norh.NoRHClient._can_handle_query(time)
    assert ans2 is False


@given(time_attr())
def test_query(time):
    qr1 = LCClient.query(time, Instrument('norh'))
    assert isinstance(qr1, QueryResponse)
    assert qr1.time_range().start == time.start
    assert qr1.time_range().end == time.end


@pytest.mark.online
@pytest.mark.parametrize("time,instrument", [
    (Time('2012/11/27', '2012/11/27'), Instrument('norh')),
    (Time('2012/10/4', '2012/10/6'), Instrument('norh')),
])
def test_get(time, instrument):
    qr1 = LCClient.query(time, instrument)
    res = LCClient.get(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@pytest.mark.online
@pytest.mark.parametrize(
    "time, instrument",
    [(a.Time('2012/10/4', '2012/10/6'), a.Instrument('norh')),
     (a.Time('2013/10/5', '2013/10/7'), a.Instrument('norh'))])
def test_fido(time, instrument):
    qr = Fido.search(time, instrument)
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
