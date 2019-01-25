"""
Provide a set of Hypothesis Strategies for various Fido related tests.
"""
import datetime

import hypothesis.strategies as st
from hypothesis import assume
from hypothesis.strategies import one_of, datetimes, sampled_from

import astropy.time
import astropy.units as u
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
    abs max is about 100 weeks + 100 days + 100 hours + a bit
    """
    keys = st.sampled_from(['days', 'seconds', 'microseconds', 'milliseconds',
                            'minutes', 'hours', 'weeks'])
    values = st.floats(min_value=1, max_value=100)
    delta = datetime.timedelta(**draw(st.dictionaries(keys, values)))
    # We don't want a 0 timedelta
    delta = astropy.time.TimeDelta(delta, format='datetime')
    assume(delta.sec > 0)
    return delta


def offline_instruments():
    """
    Returns a strategy for any instrument that does not need the internet to do
    a query
    """
    offline_instr = ['lyra', 'noaa-indices', 'noaa-predict', 'goes']
    offline_instr = st.builds(a.Instrument, st.sampled_from(offline_instr))

    return st.one_of(offline_instr)


def online_instruments():
    """
    Returns a strategy for any instrument that does not need the internet to do
    a query
    """
    online_instr = ['rhessi']
    online_instr = st.builds(a.Instrument, st.sampled_from(online_instr))

    return online_instr


@st.composite
def time_attr(draw, time=Times(
    max_value=datetime.datetime(datetime.datetime.utcnow().year, 1, 1, 0, 0),
    min_value=datetime.datetime(1900, 1, 1, 0, 0)
    ),
              delta=TimeDelta()):
    """
    Create an a.Time where it's always positive and doesn't have a massive time
    delta.
    """
    t1 = draw(time)
    t2 = t1 + draw(delta)
    # We can't download data from the future...
    assume(t2 < Time.now())

    return a.Time(t1, t2)


@st.composite
def goes_time(draw, time=Times(
    max_value=datetime.datetime(datetime.datetime.utcnow().year, 1, 1, 0, 0),
    min_value=datetime.datetime(1981, 1, 1, 0, 0)),
              delta=TimeDelta()):
    """
    Create an a.Time where it's always positive and doesn't have a massive time
    delta.
    """
    t1 = draw(time)
    t2 = t1 + draw(delta)
    # We can't download data from the future.
    assume(t2 < Time.now())

    tr = TimeRange(t1, t2)
    # There is no GOES data for this date.
    # We seem to pick up an extra micro or millisecond when we create t1
    # This seems to be test dependant.
    # `test_fido_indexing` adds 1 millisecond.
    # `test_offline_fido` adds 1 microsecond.
    # TODO: Track this issue down.
    # To avoid this we change the start of the timerange
    tr._t1 -= 1*u.s
    assume(Time('1983-05-01') not in tr)
    assume(Time('1983-05-01') + draw(delta) not in tr)

    return a.Time(tr)


def range_time(min_date, max_date=Time.now()):
    time = Times(
        min_value=datetime.datetime(1960, 1, 1, 0, 0),
        max_value=datetime.datetime(datetime.MAXYEAR, 1, 1, 0, 0),
    )
    time = time.filter(lambda x: min_date < x < max_date)
    return time_attr(time=time)
