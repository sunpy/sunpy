from __future__ import absolute_import

import numpy as np
import copy
import datetime

from sunpy.time import TimeRange
from sunpy.instr import goes
import sunpy.lightcurve as lc
from pandas.util.testing import assert_frame_equal

def test_goes_event_list():
    trange=TimeRange('2011-06-07 00:00','2011-06-08 00:00')
    result=goes.get_goes_event_list(trange,goes_class_filter='M1')
    assert type(result) == list
    assert type(result[0]) == dict
    assert type(result[0]['event_date'] == str)
    assert type(result[0]['goes_location'] == tuple)
    assert type(result[0]['peak_time'] == datetime)
    assert type(result[0]['start_time'] == datetime)
    assert type(result[0]['end_time'] == datetime)
    assert type(result[0]['goes_class'] == str)
    assert type(result[0]['noaa_active_region'] == int)
    assert result[0]['event_date'] == '2011-06-07'
    assert result[0]['goes_location'] == (54, -21)
    assert result[0]['start_time'] == datetime.datetime(2011,6,7,6,16)
    assert result[0]['peak_time'] == datetime.datetime(2011,6,7,6,41)
    assert result[0]['end_time'] == datetime.datetime(2011,6,7,6,59)
    assert result[0]['goes_class'] == 'M2.5'
    assert result[0]['noaa_active_region'] == 11226

def test_temp_em():
    # Create GOESLightcurve object, then create new one with
    # temperature & EM using with temp_em().
    goeslc = lc.GOESLightCurve.create("2014-01-01 00:00", "2014-01-01 01:00")
    goeslc_new = goes.temp_em(goeslc)
    # Find temperature and EM manually with goes_chianti_tem()
    t, em = goes.goes_chianti_tem(np.array(goeslc.data.xrsb),
                                  np.array(goeslc.data.xrsa),
                                  satellite = \
                                  int(goeslc.meta["TELESCOP"].split()[1]),
                                  date="2014-01-01")
    # Check that temperature and EM arrays from goes_chianti_tem()
    # are same as those in new GOESLightcurve object.
    assert goeslc_new.data.temperature.all() == t.all()
    assert goeslc_new.data.em.all() == em.all()
    # Check rest of data frame of new GOESLightCurve object is same
    # as that in original object.
    goeslc_revert = copy.deepcopy(goeslc_new)
    del goeslc_revert.data["temperature"]
    del goeslc_revert.data["em"]
    assert_frame_equal(goeslc_revert.data, goeslc.data)
    
def test_goes_chianti_tem():
    longflux = np.array([7e-6])
    shortflux = np.array([7e-7])
    date = "2014-04-16"
    # test case 1: satellite > 7, abundances = coronal
    temp1, em1 = goes.goes_chianti_tem(longflux, shortflux, satellite=15,
                                     date=date)
    assert np.around(temp1[0], decimals=2) == 11.28
    assert em1[0] < 4.79e+48 and em1[0] > 4.78e+48

    # test case 2: satellite > 7, abundances = photospheric
    temp2, em2 = goes.goes_chianti_tem(longflux, shortflux, satellite=15,
                                       date=date,
                                       abundances="photospheric")
        assert temp1[0] < 11.25 and temp1[0] > 11.24
    assert em2[0] < 1.12e+49 and em2[0] > 1.11e+49
    
    # test case 3: satellite < 8 and != 6, abundances = coronal
    temp3, em3 = goes.goes_chianti_tem(longflux, shortflux, satellite=5,
                                       date=date,
                                       abundances="coronal")
    assert temp3[0] < 11.43 and temp3[0] > 11.42
    assert em3[0] < 3.85e+48 and em3[0] > 3.84e+48

    # test case 4: satellite < 8 and != 6, abundances = photospheric
    temp4, em4 = goes.goes_chianti_tem(longflux, shortflux, satellite=5,
                                       date=date,
                                       abundances="photospheric")
    assert temp4[0] < 10.42 and temp4[0] > 10.41
    assert em4[0] < 8.81e+48 and em4[0] > 8.80e+48

    # test case 5: satellite = 6, date < 1983-06-28, abundances = coronal
    temp5, em5 = goes.goes_chianti_tem(longflux, shortflux, satellite=6,
                                       date="1983-06-27",
                                       abundances="coronal")
    assert temp5[0] < 12.30 and temp5[0] > 12.29
    assert em5[0] < 3.13e+48 and em5[0] > 3.12e+48
    
    # test case 6: satellite = 6, date < 1983-06-28, abundances = photospheric
    temp6, em6 = goes.goes_chianti_tem(longflux, shortflux, satellite=6,
                                       date="1983-06-27",
                                       abundances="photospheric")
    assert temp6[0] < 11.44 and temp6[0] > 11.43
    assert em6[0] < 6.74e+48 and em6[0] > 6.73e+48

    # test case 7: satellite = 6, date > 1983-06-28, abundances = coronal
    temp7, em7 = goes.goes_chianti_tem(longflux, shortflux, satellite=6,
                                       date=date,
                                       abundances="coronal")
    assert temp7[0] < 11.34 and temp7[0] > 11.33
    assert em7[0] < 4.08e+48 and em7[0] > 4.07e+48
    
    # test case 8: satellite = 6, date > 1983-06-28, abundances = photospheric
    temp8, em8 = goes.goes_chianti_tem(longflux, shortflux, satellite=6,
                                       date=date,
                                       abundances="photospheric")
    assert temp8[0] < 10.36 and temp8[0] > 10.35
    assert em8[0] < 9.94e+48 and em8[0] > 9.93e+48
    
def test_goes_get_chianti_temp():
    fluxratio = np.array([0.1])
    temp = goes.goes_get_chianti_temp(fluxratio, satellite=15)
    assert temp[0] <= 12.28 and temp >= 12.27

def test_goes_get_chianti_em():
    longflux = np.array([7e-6])
    temp = np.array([11.28])
    em = goes.goes_get_chianti_em(longflux, temp, satellite=15)
    assert em[0] <= 3.36e+48 and em[0] >= 3.35e+48
