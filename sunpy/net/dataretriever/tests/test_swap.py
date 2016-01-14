import datetime
import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument,Source,Level
from sunpy.net.dataretriever.client import QueryResponse

import sunpy.net.dataretriever.sources.swap as swap

LCClient = swap.SWAPClient()

@pytest.mark.parametrize("timerange,url_start,url_end",
                         [(TimeRange('2015/12/30 00:00:00','2015/12/30 23:59:59'),
                           'http://proba2.oma.be/swap/data/bsd/2015/12/30/swap_lv1_20151230_000044.fits',
                           'http://proba2.oma.be/swap/data/bsd/2015/12/30/swap_lv1_20151230_235935.fits')])
def test_get_url_for_time_range(timerange,url_start,url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_can_handle_query():
    ans1 = swap.SWAPClient._can_handle_query(Time('2015/12/30','2015/12/31'),Instrument('swap'))
    assert ans1 == True
    ans2 = swap.SWAPClient._can_hanlde_query(Time('2015/12/30','2015/12/31'))
    assert ans2 == False

def test_query():
    qr1 = LCClient.query(Time('2015-12-30 00:00:00','2015-12-30 00:05:00'))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 658

@pytest.mark.online
@pytest.mark.paramterize("time,instrument",
                         [(Time('2015/12/30 00:00:00','2015/12/30 00:05:00'),Instrument('swap'))])
def test_get(time,instrument):
    qr1 = LCClient.query(time,instrument)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)
                         
##urls = LCClient._get_url_for_timerange(TimeRange('2015/12/30 00:00:00','2015/12/30 23:59:59'))
##print(len(urls))
##urls.sort()
##print(urls[0])
##print(urls[-1])
