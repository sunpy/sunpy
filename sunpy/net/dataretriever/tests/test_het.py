
# This module was developed with funding from 
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

import pytest

from astropy import units as u

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.stereo as stereo

LCClient = stereo.HETClient()

@pytest.mark.parametrize("timerange, stereo_spacecraft, duration_of_average, url_start,url_end",
[(TimeRange('2008-03-01','2010-06-01'),'ahead', 15*u.min,
'http://www.srl.caltech.edu/STEREO/DATA/HET/Ahead/15minute/AeH08Apr.15m',
'http://www.srl.caltech.edu/STEREO/DATA/HET/Ahead/15minute/AeH10May.15m'),
(TimeRange('2011/06/01', '2012/09/02'),'behind', 12*u.h,
'http://www.srl.caltech.edu/STEREO/DATA/HET/Behind/12hour/BeH11Aug.12h',
'http://www.srl.caltech.edu/STEREO/DATA/HET/Behind/12hour/BeH12Sep.12h')
])

def test_get_url_for_time_range(timerange, stereo_spacecraft, duration_of_average, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange, stereo_spacecraft = stereo_spacecraft, 
                                                      duration_of_average = duration_of_average)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_can_handle_query():
    ans1 = stereo.HETClient._can_handle_query(Time(TimeRange('2007-03-01','2010-07-02')), Instrument('stereo/het'),
                                             stereo_spacecraft = 'ahead', duration_of_average = 12*u.h)
    assert ans1 == True
    ans2 = stereo.HETClient._can_handle_query(Time(TimeRange('1988/7/7', '2008/7/7')))
    assert ans2 == False
    ans3 = stereo.HETClient._can_handle_query(Time('2012/8/9', '2012/8/10'), Instrument('eve'))
    assert ans3 == False

def test_query():
    qr1 = LCClient.query(Time(TimeRange('2012/8/9', '2012/8/10')), Instrument('stereo/het'),
                                         stereo_spacecraft = 'ahead', duration_of_average = 1*u.min)
    assert isinstance(qr1,QueryResponse)
    assert len(qr1) == 1
    assert qr1.time_range()[0] == '2012/08/09'
    assert qr1.time_range()[1] == '2012/08/10'


@pytest.mark.online
@pytest.mark.parametrize("time, instrument, stereo_spacecraft, duration_of_average",
[(Time(TimeRange('2012/11/27', '2012/11/27')), Instrument('stereo/het'), stereo_spacecraft ='ahead', duration_of_average = 1*u.d),
 (Time(TimeRange('2012/10/4', '2012/10/6')), Instrument('stereo/het'), stereo_spacecraft ='behind', duration_of_average = 1*u.h),
])
def test_get(time,instrument,stereo_spacecraft,duration_of_average):
    qr1 = LCClient.query(time,instrument, stereo_spacecraft = stereo_spacecraft, duration_of_average = duration_of_average)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)
