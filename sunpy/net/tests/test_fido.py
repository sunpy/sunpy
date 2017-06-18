import os
import copy
import tempfile

import pytest
import hypothesis.strategies as st
from hypothesis import given

import astropy.units as u

from sunpy.net import attr
from sunpy.net import Fido, attrs as a
from sunpy.net.fido_factory import DownloadResponse, UnifiedResponse
from sunpy.net.dataretriever.client import CLIENTS, QueryResponse
from sunpy.util.datatype_factory_base import NoMatchError, MultipleMatchError
from sunpy.time import TimeRange, parse_time
from sunpy import config

from .strategies import (online_instruments, offline_instruments,
                         time_attr, rhessi_time, goes_time)

TIMEFORMAT = config.get("general", "time_format")


@st.composite
def offline_query(draw, instrument=offline_instruments()):
    """
    Strategy for any valid offline query
    """
    query = draw(instrument)
    # If we have AttrAnd then we don't have GOES
    if isinstance(query, a.Instrument) and query.value == 'norh':
        query &= a.Wavelength(17*u.GHz)
    if isinstance(query, a.Instrument) and query.value == 'goes':
        query &= draw(goes_time())
    else:
        query = attr.and_(query, draw(time_attr()))
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

    for block in unifiedresp.responses:
        for res in block:
            assert res.time.start in query_tr
            assert query_instr.lower() == res.instrument.lower()


@pytest.mark.online
def test_save_path():
    with tempfile.TemporaryDirectory() as target_dir:
        qr = Fido.search(a.Instrument('EVE'), a.Time("2016/10/01", "2016/10/02"), a.Level(0))
        files = Fido.fetch(qr, path=os.path.join(target_dir, "{instrument}"+os.path.sep+"{level}"))
        for f in files:
            assert target_dir in f
            assert "eve{}0".format(os.path.sep) in f


"""
Factory Tests
"""


@pytest.mark.online
def test_unified_response():
    start = parse_time("2012/1/1")
    end = parse_time("2012/1/2")
    qr = Fido.search(a.Instrument('EVE'), a.Level(0), a.Time(start, end))
    assert qr.file_num == 2
    strings = ['eve', 'SDO', start.strftime(TIMEFORMAT), end.strftime(TIMEFORMAT)]
    assert all(s in qr._repr_html_() for s in strings)


def test_no_time_error():
    query = (a.Instrument('EVE'), a.Level(0))
    with pytest.raises(ValueError) as excinfo:
        Fido.search(*query)
    assert all(str(a) in str(excinfo.value) for a in query)

    query1 = (a.Instrument('EVE') & a.Level(0))
    query2 = (a.Time("2012/1/1", "2012/1/2") & a.Instrument("AIA"))
    with pytest.raises(ValueError) as excinfo:
        Fido.search(query1 | query2)
    assert all(str(a) in str(excinfo.value) for a in query1.attrs)
    assert all(str(a) not in str(excinfo.value) for a in query2.attrs)


def test_no_match():
    with pytest.raises(NoMatchError):
        Fido.search(a.Time("2016/10/01", "2016/10/02"), a.jsoc.Series("bob"),
                    a.vso.Sample(10*u.s))


def test_call_error():
    with pytest.raises(TypeError) as excinfo:
        Fido()
    # Explicitly test all this error message as it's a copy of the one in
    # Python core.
    assert "'UnifiedDownloaderFactory' object is not callable" in str(excinfo.value)


def test_multiple_match():
    """
    Using the builtin clients a multiple match is not possible so we create a
    dummy class.
    """
    new_registry = copy.deepcopy(Fido.registry)
    Fido.registry = new_registry

    class DummyClient():
        @classmethod
        def _can_handle_query(cls, *query):
            return True
    Fido.registry.update({DummyClient: DummyClient._can_handle_query})

    with pytest.raises(MultipleMatchError):
        Fido.search(a.Time("2016/10/1", "2016/10/2"), a.Instrument('lyra'))

    Fido.registry = CLIENTS


@pytest.mark.online
def test_no_wait_fetch():
        qr = Fido.search(a.Instrument('EVE'),
                         a.Time("2016/10/01", "2016/10/02"),
                         a.Level(0))
        res = Fido.fetch(qr, wait=False)
        assert isinstance(res, DownloadResponse)
        assert isinstance(res.wait(), list)


"""
UnifiedResponse Tests

Use LYRA here because it does not use the internet to return results.
"""


def test_unifiedresponse_slicing():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument("lyra"))
    assert isinstance(results[0:2], UnifiedResponse)
    assert isinstance(results[0], UnifiedResponse)


def test_responses():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument("lyra"))

    for i, resp in enumerate(results.responses):
        assert isinstance(resp, QueryResponse)

    assert i + 1 == len(results)


def test_repr():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument("lyra"))

    rep = repr(results)
    rep = rep.split('\n')
    # 6 header lines, the results table and a blank line at the end
    assert len(rep) == 6 + len(list(results.responses)[0]) + 1
