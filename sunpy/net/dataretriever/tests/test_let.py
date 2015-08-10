
# This module was developed with funding from 
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

import pytest

from astropy import units as u

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.stereo as stereo

LCClient = stereo.LETClient()

@pytest.mark.parametrize("timerange, specie, duration_of_average, stereo_spacecraft, type_of_data, url_start, url_end",
[(TimeRange('2008-03-01','2010-07-02'),'Al', 10*u.min, 'ahead', 'summed',
'http://www.srl.caltech.edu/STEREO/DATA/Level1/Public/ahead/10Minute/2008/Summed/Al/Al_summed_ahead_2008_03_10min_level1_11.txt',
'http://www.srl.caltech.edu/STEREO/DATA/Level1/Public/ahead/10Minute/2010/Summed/Al/Al_summed_ahead_2010_07_10min_level1_11.txt'),
(TimeRange('2008/03/01', '2010/07/02'),'CNO_hi', 10*u.min, 'behind', 'sectored',
'http://www.srl.caltech.edu/STEREO/DATA/Level1/Public/behind/10Minute/2008/Sectored/CNO_hi/CNO_hi_sectored_behind_2008_03_10min_level1_11.txt',
'http://www.srl.caltech.edu/STEREO/DATA/Level1/Public/behind/10Minute/2010/Sectored/CNO_hi/CNO_hi_sectored_behind_2010_07_10min_level1_11.txt')
])

def test_get_url_for_time_range(timerange, specie, duration_of_average, stereo_spacecraft, type_of_data, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange, specie = specie, duration_of_average =duration_of_average, 
                                                        stereo_spacecraft =stereo_spacecraft, type_of_data= type_of_data)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

def test_can_handle_query():
    ans1 = stereo.LETClient._can_handle_query(Time(TimeRange('2008-03-01','2010-07-02')), Instrument('stereo/let'), specie = 'Al', 
                                                    duration_of_average = 10*u.min, stereo_spacecraft ='ahead', type_of_data ='summed')
    assert ans1 == True
    ans1 = stereo.LETClient._can_handle_query(Time(TimeRange('1998-03-01','2003-07-02')), Instrument('stereo/let'), specie = 'CNO_hi', 
                                                    duration_of_average = 10*u.min, stereo_spacecraft ='behind', type_of_data ='sectored')
    assert ans1 == True
    ans2 = stereo.LETClient._can_handle_query(Time(TimeRange('2012/7/7', '2012/7/7')))
    assert ans2 == False
    ans3 = stereo.LETClient._can_handle_query(Time(TimeRange('2012/8/9', '2012/8/10')), Instrument('eve'))
    assert ans3 == False

def test_query():
    qr1 = LCClient.query(Time(TimeRange('2012/8/9', '2012/10/11')), Instrument('stereo/let'), specie = 'Al', 
                                                    duration_of_average = 10*u.min, stereo_spacecraft ='ahead', type_of_data ='summed')
    assert isinstance(qr1,QueryResponse)
    assert len(qr1) == 2
    assert qr1.time_range()[0] == '2012/08/09'
    assert qr1.time_range()[1] == '2012/10/11'


@pytest.mark.online
@pytest.mark.parametrize("time, instrument, specie, duration_of_average, stereo_spacecraft, type_of_data",
[(Time(TimeRange('2012/01/27', '2012/04/27')), Instrument('stereo/let'), 'Al', 10*u.min, 'ahead', 'summed'),
 (Time(TimeRange('2012/10/4', '2012/12/6')), Instrument('stereo/let'),'CNO_hi', 10*u.min, 'behind', 'sectored'),
])
def test_get(time,instrument,specie,duration_of_average,stereo_spacecraft,type_of_data):
    qr1 = LCClient.query(time,instrument,specie =specie, duration_of_average = duration_of_average, 
                                        stereo_spacecraft = stereo_spacecraft, type_of_data = type_of_data)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)
