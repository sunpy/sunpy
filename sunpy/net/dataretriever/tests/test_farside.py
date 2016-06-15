#This module was developed with funding provided by
#the Google Summer of Code 2016.
import datetime
import pytest

from astropy import units as u
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.downloader_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

import sunpy.net.dataretriever.sources.gong as gong

FClient = gong.FARSIDEClient()

@pytest.mark.online
@pytest.mark.parametrize("timerange, url_start, url_end",
                         [(TimeRange('2014/2/1', '2014/2/5'),
                           'http://farside.nso.edu/oQR/fqo/201402/mrfqo140201/mrfqo140201t0000.fits',
                           'http://farside.nso.edu/oQR/fqo/201402/mrfqo140204/mrfqo140204t1200.fits')])
def test_get_url_for_timerange(timerange, url_start, url_end):
    urls = FClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_can_handle_query():
    assert FClient._can_handle_query(Time('2015/12/28', '2015/12/30'), Instrument('farside'))
    assert FClient._can_handle_query(Time('2015/12/28', '2015/12/30'), Instrument('Farside'))
    assert not FClient._can_handle_query(Time('2015/12/28', '2015/12/30'))
    assert not FClient._can_handle_query(Time('2015/12/28', '2015/12/30'), Instrument('bbso'))

@pytest.mark.online
def test_query():
    qr = FClient.query(Time('2016/1/1', '2016/1/5'), instrument = 'farside')
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 8
    assert qr.time_range()[0] == '2016/01/01'
    assert qr.time_range()[1] == '2016/01/05'

#Downloads 7 fits files each of size
#160KB. Total size ~ 1.2MB
@pytest.mark.online
@pytest.mark.parametrize("time, instrument",
                         [(Time('2016/1/1 00:00:00', '2016/1/4 10:00:00'),
                           Instrument('farside'))])
def test_get(time, instrument):
    qr = FClient.query(time, instrument)
    res = FClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)

#Downloads 5 fits files each of size 160KB.
#Total size ~ 800KB.
@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2016/5/18', '2016/5/20'), a.Instrument('farside'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
