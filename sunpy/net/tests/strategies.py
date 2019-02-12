"""
Provide a set of Hypothesis Strategies for various Fido related tests.
"""
import hypothesis.strategies as st
from hypothesis import assume
from hypothesis.strategies import datetimes

import datetime
from sunpy.net import attrs as a
from sunpy.time import TimeRange


@st.composite
def timedelta(draw):
    """
    Timedelta strategy that limits the maximum timedelta to being positive and
    abs max is about 10 weeks + 10 days + 10 hours + 10 minutes + a bit
    """
    keys = st.sampled_from(['weeks', 'days', 'hours', 'minutes', 'seconds'])
    values = st.floats(min_value=1, max_value=10)
    delta = datetime.timedelta(**draw(st.dictionaries(keys, values)))
    # We don't want a 0 timedelta
    assume(delta.total_seconds() > 0)
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
def time_attr(draw, time=datetimes(
              max_value=datetime.datetime(datetime.datetime.utcnow().year, 1, 1, 0, 0),
              min_value=datetime.datetime(1981, 1, 1, 0, 0)),
              delta=timedelta()):
    """
    Create an a.Time where it's always positive and doesn't have a massive time
    delta.
    """
    t1 = draw(time)
    t2 = t1 + draw(delta)
    # We can't download data from the future...
    assume(t2 < datetime.datetime.utcnow())

    return a.Time(t1, t2)


@st.composite
def goes_time(draw, time=datetimes(
              max_value=datetime.datetime(datetime.datetime.utcnow().year, 1, 1, 0, 0),
              min_value=datetime.datetime(1981, 1, 1, 0, 0)),
              delta=timedelta()):
    """
    Create an a.Time where it's always positive and doesn't have a massive time
    delta.
    """
    t1 = draw(time)
    delta = draw(delta)
    t2 = t1 + delta
    # We can't download data from the future.
    assume(t2 < datetime.datetime.utcnow())

    # There is no GOES data for 1983, 5, 1.
    # This checks if the date is not in the range
    assume(not (t1 <= datetime.datetime(1983, 5, 1) <= t2))
    assume(not (t1 <= (datetime.datetime(1983, 5, 1) + delta) <= t2))
    # This checks if the range start and stops on that day.
    assume((datetime.date(1983, 5, 1) - t1.date()) != datetime.timedelta(0))
    assume((datetime.date(1983, 5, 1) - t2.date()) != datetime.timedelta(0))

    tr = TimeRange(t1, t2)

    return a.Time(tr)


def range_time(min_date, max_date=datetime.datetime.utcnow()):
    time = datetimes(
        min_value=datetime.datetime(1981, 1, 1, 0, 0),
        max_value=datetime.datetime(datetime.datetime.utcnow().year, 1, 1, 0, 0),
    )
    time = time.filter(lambda x: min_date < x < max_date)
    return time_attr(time=time)
