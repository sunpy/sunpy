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

import sunpy.net.dataretriever.sources.ace as ace

SWEPAMClient = ace.SWEPAMClient()
EPAMClient = ace.EPAMClient()
MAGClient = ace.MAGClient()
SClient = ace.SISClient()

@pytest.mark.online
@pytest.mark.parametrize("timerange, url_start, url_end",
                         [(TimeRange('2015/12/27', '2015/12/30'),
                           'ftp://ftp.swpc.noaa.gov/pub/lists/ace/20151227_ace_swepam_1m.txt',
                           'ftp://ftp.swpc.noaa.gov/pub/lists/ace/20151230_ace_swepam_1m.txt')])
def test_get_url_for_timerange(timerange, url_start, url_end):
    urls = SWEPAMClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_can_handle_query():
    assert ace.SWEPAMClient._can_handle_query(Time('2015/12/30', '2015/12/31'), Instrument('swepam')) is True
    assert ace.EPAMClient._can_handle_query(Time('2015/12/30', '2015/12/31'), Instrument('epam')) is True
    assert ace.MAGClient._can_handle_query(Time('2015/12/30', '2015/12/31'), Instrument('mag')) is True
    assert ace.SISClient._can_handle_query(Time('2015/12/30', '2015/12/31'), Instrument('sis')) is True
    assert ace.SWEPAMClient._can_handle_query(Time('2015/12/30', '2015/12/31'), Instrument('swap')) is False
    assert ace.EPAMClient._can_handle_query(Time('2015/12/30', '2015/12/31')) is False

@pytest.mark.online
def test_query():
    qr = SWEPAMClient.query(Time('2015/12/27', '2015/12/30'), Instrument = 'swepam')
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 4
    assert qr.time_range()[0] == '2015/12/27'
    assert qr.time_range()[1] == '2015/12/30'

@pytest.mark.online
@pytest.mark.parametrize("time, instrument",
[(Time('2015/12/27', '2015/12/30'), Instrument('swepam'))])
def test_get(time, instrument):
    qr = SWEPAMClient.query(time, instrument)
    res = SWEPAMClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)

@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2015/12/27', '2015/12/30'), a.Instrument('epam'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
