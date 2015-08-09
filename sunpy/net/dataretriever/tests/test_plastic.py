
# This module was developed with funding from 
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

import pytest

from astropy import units as u

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.stereo as stereo

LCClient = stereo.PLASTICClient()

@pytest.mark.parametrize("timerange, stereo_spacecraft, duration_of_average, url_start,url_end",
[(TimeRange('2009-03-01','2009-04-02'),'ahead', 10*u.min,
'http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/10min/A/2009/STA_L2_PLA_1DMax_10min_20090301_060_V09.txt',
'http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/10min/A/2009/STA_L2_PLA_1DMax_10min_20090402_092_V09.txt'),
(TimeRange('2009/03/01', '2009/04/02'),'behind', 1*u.h,
'http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/1hr/B/2009/STB_L2_PLA_1DMax_1hr_20090301_060_V09.txt',
'http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/1hr/B/2009/STB_L2_PLA_1DMax_1hr_20090402_092_V09.txt')
])

def test_get_url_for_time_range(timerange, stereo_spacecraft, duration_of_average, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange, stereo_spacecraft = stereo_spacecraft, duration_of_average = duration_of_average)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_can_handle_query():
    ans1 = stereo.PLASTICClient._can_handle_query(Time(TimeRange('2009-03-01','2009-06-02')), Instrument('stereo/plastic'), 
                                                                stereo_spacecraft = 'ahead', duration_of_average = 10*u.min)
    assert ans1 == True
    ans1 = stereo.PLASTICClient._can_handle_query(Time(TimeRange('2010-03-01','2010-07-02')), Instrument('stereo/plastic'),  
                                                                stereo_spacecraft = 'behind', duration_of_average = 1*u.h)
    assert ans1 == True
    ans2 = stereo.PLASTICClient._can_handle_query(Time(TimeRange('2012/7/7', '2012/7/7')))
    assert ans2 == False
    ans3 = stereo.PLASTICClient._can_handle_query(Time(TimeRange('2012/8/9', '2012/8/10')), Instrument('eve'))
    assert ans3 == False

def test_query():
    qr1 = LCClient.query(Time('2012/8/9', '2012/8/10'), Instrument('stereo/het'), stereo_spacecraft = 'ahead', duration_of_average = 10*u.min)
    assert isinstance(qr1,QueryResponse)
    assert len(qr1) == 2
    assert qr1.time_range()[0] == '2012/08/09'
    assert qr1.time_range()[1] == '2012/08/10'


@pytest.mark.online
@pytest.mark.parametrize("time, instrument, stereo_spacecraft, duration_of_average",
[(Time(TimeRange('2012/11/27', '2012/11/27')), Instrument('stereo/het'), 'ahead', 10*u.min),
 (Time(TimeRange('2012/10/4', '2012/10/6')), Instrument('stereo/het'),'behind', 1*u.h),
])
def test_get(time,instrument, stereo_spacecraft, duration_of_average):
    qr1 = LCClient.query(time,instrument,stereo_spacecraft = stereo_spacecraft, duration_of_average = duration_of_average)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)
