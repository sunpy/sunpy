import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.goes as goes

LCClient = goes.GOESClient()

'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits'
'http://umbra.nascom.nasa.gov/goes/fits/2008/go1020080601.fits'



@pytest.mark.parametrize("timerange,url_start,url_end",
[(TimeRange('1995/06/03', '1995/06/04'),
'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits',
'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits'),
(TimeRange('2008/06/01', '2008/06/02'),
'http://umbra.nascom.nasa.gov/goes/fits/2008/go1020080601.fits',
'http://umbra.nascom.nasa.gov/goes/fits/2008/go1020080601.fits')
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

def test_can_handle_query():
    ans1 = goes.GOESClient._can_handle_query(Time('2012/8/9','2012/8/10'),Instrument('goes'))
    assert ans1 == True
    ans2 = goes.GOESClient._can_handle_query(Time('2012/7/7','2012/7/7'))
    assert ans2 == False
    ans3 = goes.GOESClient._can_handle_query(Time('2012/8/9','2012/8/10'),Instrument('eve'))
    assert ans3 == False

def test_query():
    qr1 = LCClient.query(Time('2012/8/9','2012/8/10'),Instrument('goes'))
    assert isinstance(qr1,QueryResponse)
    assert len(qr1) == 1
    assert qr1.time_range()[0] == '2012/08/09'
    assert qr1.time_range()[1] == '2012/08/10'


@pytest.mark.online
@pytest.mark.parametrize("time,instrument",
[(Time('2012/11/27','2012/11/27'),Instrument('goes')),
 (Time('2012/10/4','2012/10/6'),Instrument('goes')),
])
def test_get(time,instrument):
    qr1 = LCClient.query(time,instrument)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)

