import datetime
import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument,Source,Level
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.eve as eve

LCClient = eve.EVEClient()

@pytest.mark.parametrize("timerange,url_start,url_end",
[(TimeRange('2012/4/21','2012/4/21'),
'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120421_EVE_L0CS_DIODES_1m.txt',
'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120421_EVE_L0CS_DIODES_1m.txt'),
(TimeRange('2012/5/5','2012/5/6'),
'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120505_EVE_L0CS_DIODES_1m.txt',
'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120506_EVE_L0CS_DIODES_1m.txt'),
(TimeRange('2012/7/7','2012/7/14'),
'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120707_EVE_L0CS_DIODES_1m.txt',
'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120714_EVE_L0CS_DIODES_1m.txt')
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
    assert url == 'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2013/20130213_EVE_L0CS_DIODES_1m.txt'


def test_can_handle_query():
    ans1 = eve.EVEClient._can_handle_query(Time('2012/8/9','2012/8/10'),Instrument('eve'),Level(0))
    assert ans1 == True
    ans2 = eve.EVEClient._can_handle_query(Time('2012/7/7','2012/7/7'))
    assert ans2 == False
    ans3 = eve.EVEClient._can_handle_query(Time('2012/8/9','2012/8/10'),Instrument('eve'),Source('sdo'))
    assert ans3 ==False

def test_query():
    qr1 = LCClient.query(Time('2012/8/9','2012/8/10'),Instrument('eve'))
    assert isinstance(qr1,QueryResponse)
    assert len(qr1) == 2
    assert qr1.time_range()[0] == '2012/08/09'
    assert qr1.time_range()[1] == '2012/08/10'


@pytest.mark.online
@pytest.mark.parametrize("time,instrument",
[(Time('2012/11/27','2012/11/27'),Instrument('eve')),
 (Time('2012/10/4','2012/10/6'),Instrument('eve')),
])
def test_get(time,instrument):
    qr1 = LCClient.query(time,instrument)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)

