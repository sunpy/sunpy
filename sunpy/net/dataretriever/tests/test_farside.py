#This module was developed with funding provided by
#the Google Summer of Code 2016.
import datetime
import pytest

from astropy import units as u
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument,Physobs, Wavelength
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

trange = Time('2015/12/28', '2015/12/30')
def test_can_handle_query():
    assert FClient._can_handle_query(trange, Instrument('farside')) is True
    assert FClient._can_handle_query(trange, Instrument('Farside')) is True
    assert FClient._can_handle_query(trange) is not True
    assert FClient._can_handle_query(trange, Instrument('bbso')) is not True

@pytest.mark.online
def test_query():
    qr = FClient.query(Time('2016/1/1', '2016/1/5'), instrument = 'farside')
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 8
    assert qr.time_range()[0] == '2016/01/01'
    assert qr.time_range()[1] == '2016/01/05'
