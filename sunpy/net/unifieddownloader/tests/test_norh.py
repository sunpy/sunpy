import datetime
import pytest

import sunpy
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument,Source 
from sunpy.net.unifieddownloader.client import QueryResponse
import sunpy.net.unifieddownloader.sources.norh as norh

LCClient = norh.NoRHClient()

@pytest.mark.parametrize("timerange,url_start,url_end",
[(TimeRange('2012/4/21','2012/4/21'),
'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/04/tca120421',
'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/04/tca120421'
),
(TimeRange('2012/12/1','2012/12/2'),
'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/12/tca121201',
'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/12/tca121202'
),
(TimeRange('2012/7/7','2012/7/14'),
'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/07/tca120707',
'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/07/tca120714'
)
])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_fail_get_url_for_time_range():
    urls = LCClient._get_url_for_timerange(None)
    assert isinstance(urls, list)
    assert len(urls) == 0

def test_get_url_for_date():
    url = LCClient._get_url_for_date(datetime.date(2011,1,14))
    assert url =='ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2011/01/tca110114'


def test_can_handle_query():
    ans1 = norh.NoRHClient._can_handle_query(Time('2012/8/9','2012/8/10'),Instrument('norh'))
    assert ans1 == True
    ans2 = norh.NoRHClient._can_handle_query(Time('2012/7/7','2012/7/7'))
    assert ans2 == False

def test_query():
    qr1 = LCClient.query(Time('2012/8/9','2012/8/10'),Instrument('norh'))
    assert isinstance(qr1,QueryResponse)
    assert qr1.file_num == 2
    assert qr1.time_range()[0] == '2012/08/09'
    assert qr1.time_range()[1] == '2012/08/10'
    

@pytest.mark.online
@pytest.mark.parametrize("time,instrument",
[(Time('2012/11/27','2012/11/27'),Instrument('norh')),
 (Time('2012/10/4','2012/10/6'),Instrument('norh')),
])
def test_get(time,instrument):
    qr1 = LCClient.query(time,instrument)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == qr1.file_num

