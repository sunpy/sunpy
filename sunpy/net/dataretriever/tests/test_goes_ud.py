import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net.dataretriever.downloader_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

LCClient = goes.GOESClient()

@pytest.mark.parametrize("timerange,url_start,url_end",
[(TimeRange('1995/06/03', '1995/06/05'),
'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits',
'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950605.fits'),
(TimeRange('2008/06/02', '2008/06/04'),
'http://umbra.nascom.nasa.gov/goes/fits/2008/go1020080602.fits',
'http://umbra.nascom.nasa.gov/goes/fits/2008/go1020080604.fits')
])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_can_handle_query():
    ans1 = goes.GOESClient._can_handle_query(Time('2012/8/9', '2012/8/10'), Instrument('goes'))
    assert ans1 == True
    ans2 = goes.GOESClient._can_handle_query(Time('2012/7/7', '2012/7/7'))
    assert ans2 == False
    ans3 = goes.GOESClient._can_handle_query(Time('2012/8/9', '2012/8/10'), Instrument('eve'))
    assert ans3 == False

@pytest.mark.online
def test_query():
    qr1 = LCClient.query(Time('2012/8/9','2012/8/10'), Instrument('goes'))
    assert isinstance(qr1,QueryResponse)
    assert len(qr1) == 2
    assert qr1.time_range()[0] == '2012/08/09'
    assert qr1.time_range()[1] == '2012/08/10'


@pytest.mark.online
@pytest.mark.parametrize("time, instrument",
[(Time('2012/01/17', '2012/01/19'), Instrument('goes')),
 (Time('2012/10/4', '2012/10/6'), Instrument('goes')),
])
def test_get(time, instrument):
    qr1 = LCClient.query(time, instrument)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)

@pytest.mark.online
def test_new_logic():
    qr = LCClient.query(Time('2012/10/4','2012/10/6'), Instrument('goes'))
    res = LCClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)

@pytest.mark.online
@pytest.mark.parametrize("time, instrument",
                         [(a.Time("2012/10/4","2012/10/6"), a.Instrument("goes")),
                          (a.Time('2013/10/5', '2013/10/7'), a.Instrument("goes"))])
def test_fido(time, instrument):
    qr = Fido.search(a.Time('2012/10/4','2012/10/6'), Instrument('goes'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
