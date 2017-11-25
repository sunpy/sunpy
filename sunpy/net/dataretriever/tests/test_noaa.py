import pytest

from sunpy.time import parse_time
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.noaa as noaa
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
LCClient = noaa.NOAAIndicesClient()


@pytest.mark.parametrize(
    "timerange,url_start,url_end",
    [(TimeRange('1995/06/03', '1995/06/04'),
      'ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt',
      'ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt'),
     (TimeRange('2008/06/01', '2008/06/02'),
      'ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt',
      'ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt')])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


def test_can_handle_query():
    ans1 = noaa.NOAAIndicesClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('noaa-indices'))
    assert ans1 == True
    ans2 = noaa.NOAAIndicesClient._can_handle_query(
        Time('2012/7/7', '2012/7/7'))
    assert ans2 == False
    ans3 = noaa.NOAAIndicesClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'))
    assert ans3 == False


@pytest.mark.remote_data
def test_query():
    qr1 = LCClient.search(
        Time('2012/8/9', '2012/8/10'), Instrument('noaa-indices'))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 1
    assert qr1.time_range().start == parse_time('2012/08/09')
    assert qr1.time_range().end == parse_time('2012/08/10')


@pytest.mark.remote_data
@pytest.mark.parametrize("time, instrument", [
    (Time('2012/11/27', '2012/11/27'), Instrument('noaa-indices')),
    (Time('2012/10/4', '2012/10/6'), Instrument('noaa-indices')),
])
def test_get(time, instrument):
    qr1 = LCClient.search(time, instrument)
    res = LCClient.get(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "time, instrument",
    [(a.Time("2012/10/4", "2012/10/6"), a.Instrument('noaa-indices')),
     (a.Time('2013/10/5', '2013/10/7'), a.Instrument('noaa-indices'))])
def test_fido(time, instrument):
    qr = Fido.search(time, instrument)
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
