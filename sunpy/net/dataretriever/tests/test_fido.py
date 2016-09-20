import hypothesis.strategies as st
from hypothesis import given, assume
from hypothesis.extra.datetime import datetimes

import datetime
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


@st.composite
def time_attr(draw, time=datetimes(timezones=[]), delta=timedelta()):
    """
    Create an a.Time where it's always positive and doesn't have a massive time
    delta.
    """
    t1 = draw(time)
    t2 = t1 + draw(delta)

    return a.Time(t1, t2)


@st.composite
def goes_time(draw,
              time=datetimes(timezones=[], min_year=1981,
                             # Set Max year to prevent over large search space
                             max_year=datetime.datetime.utcnow().year),
              delta=timedelta()):
    time = time.filter(lambda x: x > parse_time('1983-04-30') and
                                 x < parse_time('1983-05-02'))
    t1 = draw(time)
    t2 = t1 + draw(delta)
    assume(t2 < datetime.datetime.utcnow())

    return a.Time(t1, t2)


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


@given(offline_query())
def test_offline_fido(query):
    unifiedresp = Fido.search(query)

    query_tr = None
    for at in query.attrs:
        if isinstance(at, a.Time):
            query_tr = TimeRange(at.start, at.end)
    if not query_tr:
        raise ValueError("No Time Specified")

    for block in unifiedresp:
        for res in block:
            assert res.time.start in query_tr
