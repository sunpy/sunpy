
# This module was developed with funding from 
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.soho as soho

LCClient = soho.ERNEClient()


@pytest.mark.parametrize("timerange,specie,url_start,url_end",
[(TimeRange('1998-03-01','2003-07-02'),'alpha',
'http://srl.utu.fi/erne_data/carrot/1933/cr1933a.txt',
'http://srl.utu.fi/erne_data/carrot/2004/cr2004a.txt'),
(TimeRange('2004/06/01', '2007/06/02'),'proton',
'http://srl.utu.fi/erne_data/carrot/2017/cr2017p.txt',
'http://srl.utu.fi/erne_data/carrot/2057/cr2057p.txt')
])

def test_get_url_for_time_range(timerange, specie, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange, specie = specie)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_can_handle_query():
    ans1 = soho.ERNEClient._can_handle_query(Time(TimeRange('1998-03-01','2003-07-02')), Instrument('soho/erne'), specie ='alpha')
    assert ans1 == True
    ans2 = soho.ERNEClient._can_handle_query(Time(TimeRange('2004-03-01','2005-07-02')), Instrument('soho/erne'), specie ='proton')
    assert ans2 == True
    ans3 = soho.ERNEClient._can_handle_query(Time(TimeRange('2012/8/9', '2012/8/10')), Instrument('eve'))
    assert ans3 == False

def test_query():
    qr1 = LCClient.query(Time(TimeRange('2006/8/9', '2008/8/10')), Instrument('soho/erne'), specie ='alpha')
    assert isinstance(qr1,QueryResponse)
    assert len(qr1) == 1
    assert qr1.time_range()[0] == '2006/08/09'
    assert qr1.time_range()[1] == '2008/08/10'


@pytest.mark.online
@pytest.mark.parametrize("time, instrument",
[(Time(TimeRange('1998-03-01','2003-07-02')), Instrument('soho/erne'),specie = 'alpha'),
 (Time(TimeRange('2004/06/01', '2007/06/02')), Instrument('soho/erne'),specie ='proton'),
])
def test_get(time,instrument,specie):
    qr1 = LCClient.query(time,instrument,specie)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)


