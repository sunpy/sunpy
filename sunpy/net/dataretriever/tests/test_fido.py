import pytest
import hypothesis.strategies as st
from hypothesis import given, assume
from hypothesis.extra.datetime import datetimes

import datetime
from functools import partial
from sunpy.net import Fido, attrs as a
from sunpy.time import TimeRange, parse_time


@st.composite
def timedelta(draw):
    """
    Timedelta strategy that limits the maximum timedelta to being positive and
    abs max is about 100 weeks + 100 days + 100 hours + a bit
    """
    keys = st.sampled_from(['days', 'seconds', 'microseconds', 'milliseconds',
                            'minutes', 'hours', 'weeks'])
    values = st.floats(min_value=0, max_value=100)
    return datetime.timedelta(**draw(st.dictionaries(keys, values)))


def offline_instruments():
    """
    Returns a strategy for any instrument that does not need the internet to do
    a query
    """
    offline_instr = ['lyra', 'norh', 'noaa-indices', 'noaa-predict', 'goes']
    offline_instr = st.builds(a.Instrument, st.sampled_from(offline_instr))

    eve = st.just(a.Instrument('eve') & a.Level(0))
    return st.one_of(offline_instr, eve)


def online_instruments():
    """
    Returns a strategy for any instrument that does not need the internet to do
    a query
    """
    online_instr = ['rhessi']
    online_instr = st.builds(a.Instrument, st.sampled_from(online_instr))

    return online_instr


@st.composite
def time_attr(draw, time=datetimes(timezones=[],
                                   max_year=datetime.datetime.utcnow().year,
                                   min_year=1900),
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


def goes_time():
    time = datetimes(timezones=[], max_year=datetime.datetime.utcnow().year,
                     min_year=1900)
    time = time.filter(lambda x: x > parse_time('1983-04-30') and
                                 x < parse_time('1983-05-02'))
    return time_attr(time=time)


def rhessi_time():
    time = datetimes(timezones=[], max_year=datetime.datetime.utcnow().year,
                     min_year=1900)
    time = time.filter(lambda x: x > parse_time('2002-02-01'))
    return time_attr(time=time)


@st.composite
def offline_query(draw, instrument=offline_instruments()):
    """
    Strategy for any valid offline query
    """
    query = draw(instrument)
    # If we have AttrAnd then we don't have GOES
    if isinstance(query, a.Instrument) and query.value == 'goes':
        query = query & draw(goes_time())
    else:
        query = query & draw(time_attr())
    return query


@st.composite
def online_query(draw, instrument=online_instruments(), time=time_attr()):
    query = draw(instrument)
    # If we have AttrAnd then we don't have RHESSI
    if isinstance(query, a.Instrument) and query.value == 'rhessi':
        query = query & draw(rhessi_time())
    return query


@given(offline_query())
def test_offline_fido(query):
    unifiedresp = Fido.search(query)
    check_response(query, unifiedresp)


@pytest.mark.online
@given(online_query())
def test_online_fido(query):
    unifiedresp = Fido.search(query)
    check_response(query, unifiedresp)


def check_response(query, unifiedresp):
    """
    Common test for online or offline query
    """
    query_tr = None
    query_instr = None
    for at in query.attrs:
        if isinstance(at, a.Time):
            query_tr = TimeRange(at.start, at.end)
        elif isinstance(at, a.Instrument):
            query_instr = at.value
    if not query_tr:
        raise ValueError("No Time Specified")

    for block in unifiedresp:
        for res in block:
            assert res.time.start in query_tr
            assert query_instr.lower() == res.instrument.lower()
