
# This module was developed with funding from 
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

import pytest

from astropy import units as u

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.stereo as stereo

LCClient = stereo.SEPTClient()

@pytest.mark.parametrize("timerange, specie, duration_of_average, stereo_spacecraft, sensor_pointing, url_start, url_end",
[(TimeRange('2008-03-01','2008-03-05'),'element', 10*u.min, 'ahead', 'asun',
'http://www2.physik.uni-kiel.de/stereo/data/sept/level2/ahead/10min/2008/sept_ahead_ele_asun_2008_061_10min_l2_v03.dat',
'http://www2.physik.uni-kiel.de/stereo/data/sept/level2/ahead/10min/2008/sept_ahead_ele_asun_2008_065_10min_l2_v03.dat'),
(TimeRange('2009/03/01', '2009/03/05'),'element', 1*u.d, 'behind', 'asun',
'http://www2.physik.uni-kiel.de/stereo/data/sept/level2/behind/1d/2009/sept_behind_ele_asun_2009_060_1d_l2_v03.dat',
'http://www2.physik.uni-kiel.de/stereo/data/sept/level2/behind/1d/2009/sept_behind_ele_asun_2009_064_1d_l2_v03.dat')
])

def test_get_url_for_time_range(timerange, specie, duration_of_average, stereo_spacecraft, sensor_pointing, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange, specie = specie, duration_of_average = duration_of_average, 
                                    stereo_spacecraft = stereo_spacecraft, sensor_pointing = sensor_pointing)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_can_handle_query():
    ans1 = stereo.SEPTClient._can_handle_query(Time(TimeRange('2008-03-01','2008-03-05')), Instrument('stereo/sept'), specie = 'element', 
                        duration_of_average = 10*u.min, stereo_spacecraft = 'ahead', sensor_pointing = 'asun')
    assert ans1 == True
    ans1 = stereo.SEPTClient._can_handle_query(Time(TimeRange('2009-03-01','2009-03-05')), Instrument('stereo/sept'), specie = 'element', 
                        duration_of_average = 1*u.d, stereo_spacecraft = 'behind', sensor_pointing = 'asun')
    assert ans1 == True
    ans2 = stereo.SEPTClient._can_handle_query(Time(TimeRange('2012/7/7', '2012/7/7')))
    assert ans2 == False
    ans3 = stereo.SEPTClient._can_handle_query(Time(TimeRange('2012/8/9', '2012/8/10')), Instrument('eve'))
    assert ans3 == False

def test_query():
    qr1 = LCClient.query(Time(TimeRange('2009/03/01', '2009/03/05')), Instrument('stereo/sept'), specie = 'element', 
                        duration_of_average = 10*u.min, stereo_spacecraft = 'ahead', sensor_pointing = 'asun')
    assert isinstance(qr1,QueryResponse)
    assert len(qr1) == 10
    assert qr1.time_range()[0] == '2009/03/01'
    assert qr1.time_range()[1] == '2009/03/05'


@pytest.mark.online
@pytest.mark.parametrize("time, instrument, specie, duration_of_average, stereo_spacecraft, sensor_pointing",
[(Time(TimeRange('2008/03/01', '2008/03/01')), Instrument('stereo/sept'), 'element', 10*u.min, 'ahead', 'asun'),
 (Time(TimeRange('2009/03/01', '2009/03/05')), Instrument('stereo/sept'), 'element', 1*u.d, 'behind', 'asun'),
])
def test_get(time,instrument, specie, duration_of_average, stereo_spacecraft, sensor_pointing):
    qr1 = LCClient.query(time,instrument,specie = specie, duration_of_average = duration_of_average, 
                                    stereo_spacecraft = stereo_spacecraft, sensor_pointing = sensor_pointing)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)
