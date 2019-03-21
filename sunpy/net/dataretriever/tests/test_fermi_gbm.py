import pytest

from sunpy.time import parse_time
from sunpy.time.timerange import TimeRange
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.fermi_gbm as fermi_gbm
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.tests.strategies import time_attr
from hypothesis import given

LCClient = fermi_gbm.GBMClient()


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange,url_start,url_end",
                         [(TimeRange('2011/06/07', '2011/06/09'),
                          'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/2011/06/07/'
                          'current/glg_cspec_n5_110607_v00.pha',
                          'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/2011/06/09/'
                          'current/glg_cspec_n5_110609_v00.pha'),
                          (TimeRange('2016/09/09', '2016/09/11'),
                          'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/2016/09/09/'
                          'current/glg_cspec_n5_160909_v00.pha',
                          'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'
                          '2016/09/11/current/glg_cspec_n5_160911_v00.pha')])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange, detector='n5', resolution='cspec')
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@given(time_attr())
def test_can_handle_query(time):
    ans1 = LCClient._can_handle_query(time, a.Instrument('gbm'))
    assert ans1 is True
    ans2 = LCClient._can_handle_query(time, a.Instrument('gbm'),
                                      a.Detector('n5'))
    assert ans2 is True
    ans3 = LCClient._can_handle_query(time, a.Instrument('gbm'),
                                      a.Detector('n5'), a.Resolution('ctime'))
    assert ans3 is True
    ans4 = LCClient._can_handle_query(time)
    assert ans4 is False


@pytest.mark.remote_data
@pytest.mark.parametrize("time,instrument", [
    (a.Time('2012/8/9', '2012/8/10'), a.Instrument('gbm')),
])
def test_query(time, instrument):
    qr1 = LCClient.search(time, instrument)
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 2
    assert qr1.time_range().start == time.start
    assert qr1.time_range().end == time.end


@pytest.mark.remote_data
@pytest.mark.parametrize("time,instrument", [
    (a.Time('2012/11/27', '2012/11/27'), a.Instrument('gbm')),
])
def test_get(time, instrument):
    qr1 = LCClient.search(time, instrument)
    download_list = LCClient.fetch(qr1)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    'query',
    [(a.Time('2012/10/4', '2012/10/6') & a.Instrument('gbm') & a.Detector('n5'))])
def test_fido(query):
    qr = Fido.search(query)
    client = qr.get_response(0).client
    assert isinstance(qr, UnifiedResponse)
    assert type(client) == type(LCClient)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
