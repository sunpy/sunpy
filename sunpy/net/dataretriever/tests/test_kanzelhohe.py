#This module was developed with funding provided by
#the Google Summer of Code 2016.
import datetime
import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument,Level
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.downloader_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

import sunpy.net.dataretriever.sources.kanzelhohe as kanzelhohe

LCClient = kanzelhohe.KanzelhoheClient()

@pytest.mark.online
@pytest.mark.parametrize("timerange, url_start, url_end",
                         [(TimeRange('2015/01/10 00:00:00', '2015/01/10 12:00:00'),
                           'http://cesar.kso.ac.at/halpha2k/recent/2015/kanz_halph_fr_20150110_102629.fts.gz',
                           'http://cesar.kso.ac.at/halpha2k/recent/2015/kanz_halph_fr_20150110_113524.fts.gz')])
def test_get_url_for_timerange(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    print(urls)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[1] == url_end

def test_can_handle_query():
    ans1 = kanzelhohe.KanzelhoheClient._can_handle_query(Time('2015/12/30 00:00:00','2015/12/31 00:05:00'), Instrument('kanzelhohe'))
    assert ans1 == True
    ans2 = kanzelhohe.KanzelhoheClient._can_handle_query(Time('2015/12/30 00:00:00','2015/12/31 00:05:00'), Instrument('swap'))
    assert ans2 == False
    ans3 = kanzelhohe.KanzelhoheClient._can_handle_query(Time('2015/12/30 00:00:00','2015/12/31 00:05:00'))
    assert ans3 == False

@pytest.mark.online
def test_query():
    qr = LCClient.query(Time('2015/01/10 00:00:00', '2015/01/10 12:00:00'))
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 2
    assert qr.time_range()[0] == '2015/01/10'
    assert qr.time_range()[1] == '2015/01/10'

@pytest.mark.online
@pytest.mark.parametrize("time, instrument",
                         [(Time('2015/01/10 00:00:00', '2015/01/10 12:00:00'),
                           Instrument('kanzelhohe'))])                          
def test_get(time, instrument):
    qr = LCClient.query(time, instrument)
    res = LCClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)

@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2015/01/10 00:00:00', '2015/01/10 12:00:00'), a.Instrument('kanzelhohe') )
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
