import pytest
import hypothesis.strategies as st
from hypothesis import given

from sunpy.net import Fido, attrs as a
from sunpy.time import TimeRange

from .strategies import (online_instruments, offline_instruments,
                         time_attr, rhessi_time, goes_time)


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
