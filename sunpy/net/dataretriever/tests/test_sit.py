
# This module was developed with funding from 
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

import pytest

from astropy import units as u

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.stereo as stereo

LCClient = stereo.SITClient()

@pytest.mark.parametrize("timerange, specie, stereo_spacecraft, duration_of_average, url_start,url_end",
[(TimeRange('2008-03-01','2008-07-02'),'4He', 'ahead', 10*u.min, 
'http://www.srl.caltech.edu/STEREO/DATA/SIT/ahead/10min/4He/SIT_Ahead_10min_4He_2008_03.txt',
'http://www.srl.caltech.edu/STEREO/DATA/SIT/ahead/10min/4He/SIT_Ahead_10min_4He_2008_07.txt'),
(TimeRange('2007/04/01', '2009/09/02'),'O', 'behind', 1*u.h,
'http://www.srl.caltech.edu/STEREO/DATA/SIT/behind/1hr/SIT_Behind_1hr_O_2008.txt',
'http://www.srl.caltech.edu/STEREO/DATA/SIT/behind/1hr/SIT_Behind_1hr_O_2009.txt')
])

def test_get_url_for_time_range(timerange, specie, stereo_spacecraft, duration_of_average, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange, specie = specie, stereo_spacecraft = stereo_spacecraft,
                                                duration_of_average = duration_of_average)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_can_handle_query():
    ans1 = stereo.SITClient._can_handle_query(Time(TimeRange('2008-03-01','2008-07-02')), Instrument('stereo/sit'), specie = '4He',
                                                stereo_spacecraft = 'ahead', duration_of_average = 10*u.min)
    assert ans1 == True
    ans1 = stereo.SITClient._can_handle_query(Time(TimeRange('2007-04-01','2009-09-02')), Instrument('stereo/sit'), specie = 'O',
                                                stereo_spacecraft = 'behind', duration_of_average = 1*u.h)
    assert ans1 == True
    ans2 = stereo.SITClient._can_handle_query(Time(TimeRange('2012/7/7', '2012/7/7')))
    assert ans2 == False
    ans3 = stereo.SITClient._can_handle_query(Time(TimeRange('2012/8/9', '2012/8/10')), Instrument('eve'))
    assert ans3 == False

def test_query():
    qr1 = LCClient.query(Time(TimeRange('2008/03/01', '2008/04/02')), Instrument('stereo/sit'), specie = '4He',
                                                stereo_spacecraft = 'ahead', duration_of_average = 10*u.min)
    assert isinstance(qr1,QueryResponse)
    assert len(qr1) == 2
    assert qr1.time_range()[0] == '2008/03/01'
    assert qr1.time_range()[1] == '2008/04/02'


@pytest.mark.online
@pytest.mark.parametrize("time, instrument, specie, stereo_spacecraft, duration_of_average",
[(Time(TimeRange('2008/03/01', '2008/07/02')), Instrument('stereo/sit'),'4He', 'ahead', 10*u.min),
 (Time(TimeRange('2007/04/01', '2009/09/02')), Instrument('stereo/sit'), 'O', 'behind', 1*u.h),
])
def test_get(time, instrument, specie, stereo_spacecraft, duration_of_average ):
    qr1 = LCClient.query(time,instrument,specie = specie, stereo_spacecraft = stereo_spacecraft, duration_of_average = duration_of_average)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)
