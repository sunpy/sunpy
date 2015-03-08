import datetime
import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.unifieddownloader.client import QueryResponse
import sunpy.net.unifieddownloader.sources.lyra as lyra

LCClient = lyra.LYRAClient()

@pytest.mark.parametrize("timerange,url_start,url_end",
[
(TimeRange('2012/1/7','2012/1/7'),
'http://proba2.oma.be/lyra/data/bsd/2012/01/07/lyra_20120107-000000_lev2_std.fits',
'http://proba2.oma.be/lyra/data/bsd/2012/01/07/lyra_20120107-000000_lev2_std.fits'
),
(TimeRange('2012/12/1','2012/12/2'),
'http://proba2.oma.be/lyra/data/bsd/2012/12/01/lyra_20121201-000000_lev2_std.fits',
'http://proba2.oma.be/lyra/data/bsd/2012/12/02/lyra_20121202-000000_lev2_std.fits'
),
(TimeRange('2012/4/7','2012/4/14'),
'http://proba2.oma.be/lyra/data/bsd/2012/04/07/lyra_20120407-000000_lev2_std.fits',
'http://proba2.oma.be/lyra/data/bsd/2012/04/14/lyra_20120414-000000_lev2_std.fits')
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
    url = LCClient._get_url_for_date(datetime.date(2013,2,13))
    assert url == 'http://proba2.oma.be/lyra/data/bsd/2013/02/13/lyra_20130213-000000_lev2_std.fits'

def test_can_handle_query():
    ans1 = lyra.LYRAClient._can_handle_query(Time('2011/8/9','2011/8/10'),Instrument('lyra'))
    assert ans1 == True
    ans2 = lyra.LYRAClient._can_handle_query(Time('2013/1/7','2013/1/7'))
    assert ans2 == False


def test_query():
    qr1 = LCClient.query(Time('2012/11/10','2012/11/11'),Instrument('lyra'))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 2
    assert qr1.time_range()[0] == '2012/11/10'
    assert qr1.time_range()[1] == '2012/11/11'
    

@pytest.mark.online
@pytest.mark.parametrize("time,instrument",
[(Time('2013/8/27','2013/8/27'),Instrument('lyra')),
 (Time('2013/2/4','2013/2/6'),Instrument('lyra')),
])
def test_get(time,instrument):
    qr1 = LCClient.query(time,instrument)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)

