import pytest

from unittest import mock

from sunpy.time import parse_time
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.noaa as noaa
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
LCClient = noaa.NOAAIndicesClient()


def create_mock_unified_object(start_date, end_date):
    '''
    Let us create a mock QueryResponse object
    using that, we will construct UnifiedResponse
    We will prefill some downloaded data from
    running noaa.NOAAIndicesClient().search(Time('2012/8/9', '2012/8/10'),
                                                    Instrument('noaa-indices'))
    '''
    # Create a mock QueryResponse object
    map_ = {}
    map_['Time_start'] = parse_time(start_date)
    map_['Time_end'] = parse_time(end_date)
    map_['source'] = 'sdic'
    map_['instrument'] = 'noaa-indices'
    map_['physobs'] = 'sunspot number'
    map_['provider'] = 'swpc'

    resp = QueryResponse.create(map_, [''])
    # Attach the client with the QueryResponse
    resp.client = 'noaa-indices'

    # Create a UnifiedResponse object
    uresp = UnifiedResponse(resp)
    print(type(uresp))
    return (uresp)


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
    assert ans1 is True
    ans2 = noaa.NOAAIndicesClient._can_handle_query(
        Time('2012/7/7', '2012/7/7'))
    assert ans2 is False
    ans3 = noaa.NOAAIndicesClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'))
    assert ans3 is False


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
def test_fetch(time, instrument):
    qr1 = LCClient.search(time, instrument)
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@mock.patch('sunpy.net.fido_factory.Fido.fetch',
            side_effect=create_mock_unified_object("2012/10/4", "2012/10/6"))
def test_fido(mock_search):
    qr = Fido.search(a.Time("2012/10/4", "2012/10/6"),
                     a.Instrument('noaa-indices'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert isinstance(response, QueryResponse)
    assert len(response) == qr._numfile
