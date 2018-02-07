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
    abs max is about 100 weeks + 100 days + 100 hours + a bit
    """
    keys = st.sampled_from(['days', 'seconds', 'microseconds', 'milliseconds',
                            'minutes', 'hours', 'weeks'])
    values = st.floats(min_value=1, max_value=100)
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
    min_value=datetime.datetime(1900, 1, 1, 0, 0)
    ),
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
    t2 = t1 + draw(delta)
    # We can't download data from the future.
    assume(t2 < datetime.datetime.utcnow())

    tr = TimeRange(t1, t2)
    # There is no GOES data for this date.
    assume(datetime.datetime(1983, 5, 1, 0, 0, 0) not in tr)
    assume((datetime.datetime(1983, 5, 1) + draw(delta)) not in tr)

    return a.Time(tr)


def range_time(min_date, max_date=datetime.datetime.utcnow()):
    time = datetimes(
        min_value=datetime.datetime(1960, 1, 1, 0, 0),
        max_value=datetime.datetime(datetime.MAXYEAR, 1, 1, 0, 0),
    )
    time = time.filter(lambda x: min_date < x < max_date)
    return time_attr(time=time)
