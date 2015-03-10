import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.rhessi as rhessi

LCClient = rhessi.RHESSIClient()

@pytest.mark.parametrize("timerange,url_start",
[
(TimeRange('2012/7/1','2012/7/1'),
'http://hesperia.gsfc.nasa.gov/hessidata/metadata/catalog/hsi_obssumm_20120701_050.fits'),
(TimeRange('2013/6/3','2013/6/4'),
'http://hesperia.gsfc.nasa.gov/hessidata/metadata/catalog/hsi_obssumm_20130603_042.fits'),
(TimeRange('2012/7/1','2012/7/14'),
'http://hesperia.gsfc.nasa.gov/hessidata/metadata/catalog/hsi_obssumm_20120701_050.fits')
])
def test_get_url_for_time_range(timerange, url_start):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start

def test_fail_get_url_for_time_range():
    urls = LCClient._get_url_for_timerange(None)
    assert isinstance(urls, list)
    assert len(urls) == 0

def test_can_handle_query():
    ans1 = rhessi.RHESSIClient._can_handle_query(Time('2012/8/9','2012/8/9'),Instrument('rhessi'))
    assert ans1 == True
    ans2 = rhessi.RHESSIClient._can_handle_query(Time('2013/2/7','2013/2/7'))
    assert ans2 == False

def test_query():
    qr1 = LCClient.query(Time('2011/4/9','2011/4/9'),Instrument('rhessi'))
    assert isinstance(qr1,QueryResponse)
    assert len(qr1) == 1
    assert qr1.time_range()[0] == '2011/04/09'
    assert qr1.time_range()[1] == '2011/04/09'


@pytest.mark.online
@pytest.mark.parametrize("time,instrument",
[(Time('2012/11/27','2012/11/27'),Instrument('rhessi')),
 (Time('2012/10/4','2012/10/4'),Instrument('rhessi')),
])
def test_get(time,instrument):
    qr1 = LCClient.query(time,instrument)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)

