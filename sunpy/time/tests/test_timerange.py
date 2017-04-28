# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import datetime

import pytest
import astropy.units as u

import sunpy.time
from sunpy.extern.six.moves import zip

tbegin_str = '2012/1/1'
tfin_str = '2012/1/2'
dt = u.Quantity(24 * 60 * 60, 's')

start = sunpy.time.parse_time(tbegin_str)
end = sunpy.time.parse_time(tfin_str)
delta = end - start


@pytest.mark.parametrize("inputs", [
    (tbegin_str, tfin_str),
    (tbegin_str, dt),
    (tbegin_str, datetime.timedelta(days=1))
])
def test_timerange_inputs(inputs):
    timerange = sunpy.time.TimeRange(*inputs)
    assert isinstance(timerange, sunpy.time.TimeRange)
    assert timerange.start == start
    assert timerange.end == end
    assert timerange.dt == delta


def test_timerange_invalid_range():
    """
    TimeRange can expect a range of two times. if this is not the case then
    raise an exception
    """
    lower = '2016/01/04 09:30'
    mid = '2016/06/04 09:30'
    upper = '2017/03/04 09:30'

    with pytest.raises(ValueError):
        sunpy.time.TimeRange((lower,))

    with pytest.raises(ValueError):
        sunpy.time.TimeRange((lower, mid, upper))


def test_equals():
    """
    Test the __eq__ operator
    """
    lower = '2016/01/04T09:30:00.000'
    upper = '2016/06/04T09:30:00.000'
    upper_plus_one_msec = '2016/06/04T09:30:00.001'

    tr = sunpy.time.TimeRange((lower, upper))

    # This should *always* hold true
    assert tr == tr

    # Same values, different format
    tr_diff_format = sunpy.time.TimeRange('2016-01-04T09:30:00.000', '2016-06-04T09:30:00.000')

    assert tr == tr_diff_format

    lower_dt = datetime.datetime(year=2016, month=1, day=4, hour=9, minute=30, second=0)
    upper_dt = datetime.datetime(year=2016, month=6, day=4, hour=9, minute=30, second=0)
    tr_datetime = sunpy.time.TimeRange(lower_dt, upper_dt)

    assert tr == tr_datetime

    tr_plus_one_msec = sunpy.time.TimeRange((lower, upper_plus_one_msec))

    assert (tr_plus_one_msec == tr) is False

    # Attempt using objects which are *not* TimeRanges
    assert (tr == lower_dt) is False
    assert (lower_dt == tr) is False


def test_not_equals():
    """
    Test __ne__ operator
    """
    a_st = '2016/01/04T09:30:00.000'
    a_et = '2016/06/04T09:30:00.000'

    b_st = '2017/01/04T09:30:00.000'
    b_et = '2017/06/04T09:30:00.000'

    # Same start time, different end times
    assert sunpy.time.TimeRange(a_st, a_et) != sunpy.time.TimeRange(a_st, b_et)

    # Different start times, same end times
    assert sunpy.time.TimeRange(b_st, b_et) != sunpy.time.TimeRange(a_st, b_et)

    # Different start & end times
    assert sunpy.time.TimeRange(a_st, a_et) != sunpy.time.TimeRange(b_st, b_et)

    # Different objects
    assert sunpy.time.TimeRange(a_st, a_et) != dict()
    assert list() != sunpy.time.TimeRange(a_st, a_et)


def test_get_dates():
    lower = '2016/01/04 09:30'
    lower_plus_one_day = '2016/01/05 09:30'
    upper = '2016/06/04 09:30'

    single_day = sunpy.time.TimeRange((lower, lower))

    assert single_day.get_dates() == [datetime.date(2016, 1, 4)]

    two_days = sunpy.time.TimeRange((lower, lower_plus_one_day))

    assert two_days.get_dates() == [datetime.date(2016, 1, 4), datetime.date(2016, 1, 5)]

    one_year = sunpy.time.TimeRange('2017/01/01', '2017-12-31')
    assert len(one_year.get_dates()) == 365

    leap_year = sunpy.time.TimeRange('2016/01/01', '2016-12-31')
    assert len(leap_year.get_dates()) == 366


@pytest.mark.parametrize("ainput", [
    (tbegin_str, tfin_str),
    (tbegin_str, dt),
    (tbegin_str, datetime.timedelta(days=1)),
    (sunpy.time.TimeRange(tbegin_str, tfin_str))
    ])
def test_timerange_input(ainput):
    timerange = sunpy.time.TimeRange(ainput)
    assert isinstance(timerange, sunpy.time.TimeRange)
    assert timerange.start == start
    assert timerange.end == end
    assert timerange.dt == delta


@pytest.mark.parametrize("ainput", [
    (tbegin_str, tfin_str),
    (tfin_str, -dt),
    (tfin_str, tbegin_str)
    ])
def test_start_lessthan_end(ainput):
    """Test that the start and end time for a timerange is always in the
    right order"""
    timerange = sunpy.time.TimeRange(ainput)
    t1 = timerange.start
    t2 = timerange.end
    assert t1 < t2
    assert timerange.start == start
    assert timerange.end == end


@pytest.fixture
def timerange_a():
    return sunpy.time.TimeRange(tbegin_str, tfin_str)


def test_center(timerange_a):
    assert timerange_a.center == datetime.datetime(year=2012, day=1, month=1,
                                                   hour=12)


def test_split(timerange_a):
    expect = [sunpy.time.TimeRange('2012/1/1T00:00:00', '2012/1/1T12:00:00'),
              sunpy.time.TimeRange('2012/1/1T12:00:00', '2012/1/2T00:00:00')]
    split = timerange_a.split(n=2)
    # Doing direct comparisons seem to not work
    assert all([wi.start == ex.start and wi.end == ex.end for wi, ex in zip(split, expect)])


def test_split_n_0_error(timerange_a):
    with pytest.raises(ValueError):
        timerange_a.split(n=0)


def test_input_error(timerange_a):
    with pytest.raises(ValueError):
        sunpy.time.TimeRange((tbegin_str))


def test_window(timerange_a):
    timerange = sunpy.time.TimeRange(tbegin_str, tfin_str)
    window = timerange.window(u.Quantity(12 * 60 * 60, 's'), u.Quantity(10, 's'))
    expect = [sunpy.time.TimeRange('2012/1/1T00:00:00', '2012/1/1T00:00:10'),
              sunpy.time.TimeRange('2012/1/1T12:00:00', '2012/1/1T12:00:10'),
              sunpy.time.TimeRange('2012/1/2T00:00:00', '2012/1/2T00:00:10')]
    assert isinstance(window, list)
    # Doing direct comparisons seem to not work
    assert all([wi.start == ex.start and wi.end == ex.end for wi, ex in zip(window, expect)])


def test_window_timedelta(timerange_a):
    timerange = sunpy.time.TimeRange(tbegin_str, tfin_str)
    window = timerange.window(datetime.timedelta(hours=12), datetime.timedelta(seconds=10))
    expect = [sunpy.time.TimeRange('2012/1/1T00:00:00', '2012/1/1T00:00:10'),
              sunpy.time.TimeRange('2012/1/1T12:00:00', '2012/1/1T12:00:10'),
              sunpy.time.TimeRange('2012/1/2T00:00:00', '2012/1/2T00:00:10')]
    assert isinstance(window, list)
    # Doing direct comparisons seem to not work
    assert all([wi.start == ex.start and wi.end == ex.end for wi, ex in zip(window, expect)])


def test_days(timerange_a):
    assert timerange_a.days == u.Quantity(1, 'd')


def test_start(timerange_a):
    assert timerange_a.start == start


def test_end(timerange_a):
    assert timerange_a.end == end


def test_seconds(timerange_a):
    assert timerange_a.seconds == dt


def test_minutes(timerange_a):
    assert timerange_a.minutes == u.Quantity(24 * 60, 'min')


def test_hours(timerange_a):
    assert timerange_a.hours == u.Quantity(24, 'hour')


def test_next():
    timerange = sunpy.time.TimeRange(tbegin_str, tfin_str)
    timerange.next()
    assert isinstance(timerange, sunpy.time.TimeRange)
    assert timerange.start == start + delta
    assert timerange.end == end + delta
    assert timerange.dt == delta


def test_previous():
    timerange = sunpy.time.TimeRange(tbegin_str, tfin_str)
    timerange.previous()
    assert isinstance(timerange, sunpy.time.TimeRange)
    assert timerange.start == start - delta
    assert timerange.end == end - delta
    assert timerange.dt == delta


def test_extend():
    timerange = sunpy.time.TimeRange(tbegin_str, tfin_str)
    timerange.extend(delta, delta)
    assert isinstance(timerange, sunpy.time.TimeRange)
    assert timerange.start == start + delta
    assert timerange.end == end + delta
    assert timerange.dt == delta


def test_contains(timerange_a):
    before = datetime.datetime(year=1990, month=1, day=1)
    after = datetime.datetime(year=2022, month=1, day=1)
    between = datetime.datetime(year=2014, month=5, day=4)
    timerange = sunpy.time.TimeRange('2014/05/03 12:00', '2014/05/05 21:00')
    assert between in timerange
    assert before not in timerange
    assert after not in timerange
    assert timerange.start in timerange
    assert timerange.end in timerange
    assert '2014/05/04 15:21' in timerange
    assert '1975/4/13' not in timerange
    assert '2100/1/1'not in timerange
    assert '2014/05/03 12:00' in timerange
    assert '2014/05/05 21:00' in timerange
