"""
Provide a set of Hypothesis Strategies for various Fido related tests.
"""
import datetime

import numpy as np
import hypothesis.strategies as st
from hypothesis import assume
from hypothesis.strategies import one_of, datetimes, sampled_from

import astropy.time
from astropy.time import Time

from sunpy.net import attrs as a
from sunpy.time import TimeRange

TimesLeapsecond = sampled_from((Time('2015-06-30T23:59:60'),
                                Time('2012-06-30T23:59:60')))


@st.composite
def Times(draw, max_value, min_value):
    time = one_of(datetimes(max_value=max_value, min_value=min_value),
                  TimesLeapsecond)

    time = Time(draw(time))

    return time


@st.composite
def TimeDelta(draw):
    """
    Timedelta strategy that limits the maximum timedelta to being positive and
    abs max is about 10 weeks + 10 days + 10 hours + 10 minutes + a bit
    """
    keys = st.sampled_from(['weeks', 'days', 'hours', 'minutes', 'seconds'])
    values = st.floats(min_value=1, max_value=10)
    time_dict = {'days': st.floats(min_value=1, max_value=8),
                 'hours': st.floats(min_value=1, max_value=12),
                 'minutes': st.floats(min_value=1, max_value=30),
                 'seconds': st.floats(min_value=1, max_value=60)}

    delta = datetime.timedelta(**draw(st.fixed_dictionaries(time_dict)))
    delta = astropy.time.TimeDelta(delta, format='datetime')

    # We don't want a 0 timedelta
    assume(delta.sec > 0)

    return delta


def offline_instruments():
    """
    Returns a strategy for any instrument that does not need the internet to do
    a query.
    """
    offline_instr = ['lyra', 'noaa-indices', 'noaa-predict', 'soon', 'goes']
    offline_instr = st.builds(a.Instrument, st.sampled_from(offline_instr))

    return st.one_of(offline_instr)


def online_instruments():
    """
    Returns a strategy for any instrument that does need the internet to do
    a query.
    """
    online_instr = ['eve', 'rhessi', 'norh']
    online_instr = st.builds(a.Instrument, st.sampled_from(online_instr))

    return online_instr


@st.composite
def time_attr(draw, time=Times(
              max_value=datetime.datetime(datetime.datetime.utcnow().year, 1, 1, 0, 0),
              min_value=datetime.datetime(1981, 1, 1, 0, 0)),
              delta=TimeDelta()):
    """
    Create an a.Time where it's always positive.
    """
    t1 = draw(time)
    t2 = t1 + draw(delta)
    # We can't download data from the future.
    assume(t2 < Time.now())

    return a.Time(t1, t2)


@st.composite
def goes_time(draw, time=Times(
              max_value=datetime.datetime(datetime.datetime.utcnow().year, 1, 1, 0, 0),
              min_value=datetime.datetime(1981, 1, 1, 0, 0)),
              delta=TimeDelta()):
    """
    Create an a.Time where it's always positive.
    """
    t1 = draw(time)
    delta = draw(delta)
    t2 = t1 + delta
    # We can't download data from the future.
    assume(t2 < Time.now())

    # There is no GOES data for this date.
    assume(not (t1 <= Time('1983-05-01') <= t2))
    assume(not (t1 <= (Time('1983-05-01') + delta) <= t2))
    # This checks if the range start and stops on that day.
    assume((np.abs(Time('1983-05-01') - t1)) > astropy.time.TimeDelta(0.01))
    assume((np.abs(Time('1983-05-01') - t2)) > astropy.time.TimeDelta(0.01))

    tr = TimeRange(t1, t2)

    return a.Time(tr)


def range_time(min_date, max_date=Time.now()):
    time = Times(
        min_value=datetime.datetime(1981, 1, 1, 0, 0),
        max_value=datetime.datetime(datetime.datetime.utcnow().year, 1, 1, 0, 0),
    )
    time = time.filter(lambda x: min_date < x < max_date)

    return time_attr(time=time)
