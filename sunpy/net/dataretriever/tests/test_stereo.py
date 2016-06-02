#This module was developed with funding provided by
#the Google Summer of Code 2016.
import datetime
import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument,Source,Detector
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.downloader_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
import sunpy.net.dataretriever.sources.stereo as stereo

SECCHIClient = stereo.SECCHIClient()
PClient = stereo.PLASTICClient()
IClient = stereo.IMPACTClient()
SWClient = stereo.SWAVESClient()

@pytest.mark.online
@pytest.mark.parametrize("timerange, source, detector, url_start, url_end",
                         [(TimeRange('2008/3/2 17:45:00', '2008/3/2 18:00:00'), 'ahead', 'euvi',
                          'http://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/secchi/img/euvi/20080302/20080302_174530_n7euA.fts',
                         'http://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/secchi/img/euvi/20080302/20080302_175530_n7euA.fts')])
def test_get_url_for_timerange(timerange, source, detector, url_start, url_end):
    urls = SECCHIClient._get_url_for_timerange(timerange, source = source, detector = detector)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

trange = Time('2008/3/2 17:45:00', '2008/3/2 18:00:00')
def test_can_handle_query():
    assert SECCHIClient._can_handle_query(trange, Instrument('secchi'), Source('ahead'), Detector('euvi')) is True
    assert SECCHIClient._can_handle_query(trange, Instrument('secchi'), Source('behind'), Detector('hi_1')) is True
    assert PClient._can_handle_query(trange, Instrument('plastic'), Source('ahead')) is True
    assert IClient._can_handle_query(trange, Instrument('impact'), Source('Behind')) is True
    assert SWClient._can_handle_query(trange, Instrument('swaves')) is not True
    assert SECCHIClient._can_handle_query(trange, Instrument('secchi'), Source('Ahead'), Detector('xyz')) is not True

trange = Time('2007/6/13 04:00:00', '2007/6/13 05:00:00')
@pytest.mark.online
def test_query():
    qr = SECCHIClient.query(trange, Instrument = 'secchi', source = 'behind', detector = 'euvi')
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 5
    assert qr.time_range()[0] == '2007/06/13'
    assert qr.time_range()[1] == '2007/06/13'

@pytest.mark.online
@pytest.mark.parametrize("time, instrument, source",
[(Time('2013/4/1', '2013/4/3'), Instrument("plastic"), Source("ahead"))])
def test_get(time, instrument, source):
    qr = PClient.query(time, instrument, source)
    res = PClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)


              
