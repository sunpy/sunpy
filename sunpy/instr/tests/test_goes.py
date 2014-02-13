from __future__ import absolute_import

import datetime
from sunpy.time import TimeRange
from sunpy.instr import goes

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
