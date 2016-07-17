#This module was developed with funding provided by
#the Google Summer of Code 2016.
import datetime
import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument, Source, Detector
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
                         [(TimeRange('2008/3/2 17:45:00', '2008/3/2 18:00:00'), 'STEREO_A', 'euvi',
                          'http://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/secchi/img/euvi/20080302/20080302_174530_n7euA.fts',
                         'http://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/secchi/img/euvi/20080302/20080302_175530_n7euA.fts')])
def test_get_url_for_timerange(timerange, source, detector, url_start, url_end):
    urls = SECCHIClient._get_url_for_timerange(timerange, source = source, detector = detector)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

trange = Time('2008/3/2 17:45:00', '2008/3/2 18:00:00')
@pytest.mark.parametrize("client, timerange, instrument, source, detector, expected",
                         [(SECCHIClient, trange, Instrument('secchi'), Source('STEREO_A'), Detector('euvi'), True),
                          (SECCHIClient, trange, Instrument('secchi'), Source('STEREO_B'), Detector('hi_1'), True),
                          (PClient, trange, Instrument('plastic'), Source('STEREO_A'), None, True),
                          (IClient, trange, Instrument('impact'), Source('STEREO_B'), None, True),
                          (SWClient, trange, Instrument('swaves'), None, None, False),
                          (SECCHIClient, trange, Instrument('secchi'), Source("STEREO_A"), Detector('xyz'), False)])
def test_can_handle_query(client, timerange, instrument, source, detector, expected):
    assert client._can_handle_query(timerange, instrument, source, detector) is expected

trange = Time('2007/6/13 04:00:00', '2007/6/13 05:00:00')
@pytest.mark.online
def test_query():
    qr = SECCHIClient.query(trange, Instrument = 'secchi', source = 'STEREO_B', detector = 'euvi')
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 5
    assert qr.time_range()[0] == '2007/06/13'
    assert qr.time_range()[1] == '2007/06/13'

##Downloads 2 files each of size 1.4MB
##Total size = 2.8MB
@pytest.mark.online
@pytest.mark.parametrize("time, instrument, source",
[(Time('2013/4/1', '2013/4/3'), Instrument("plastic"), Source("STEREO_A"))])
def test_get(time, instrument, source):
    qr = PClient.query(time, instrument, source)
    res = PClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)


#Downloads 3 files each of size 149KB
#Total size is 450KB.
@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2013/4/5 01:00:00', '2013/4/5 06:00:00'), a.Instrument('secchi'), a.Source('STEREO_B'), a.Detector('hi_2'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
    
