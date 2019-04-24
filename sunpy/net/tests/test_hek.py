# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

#pylint: disable=W0613

import pytest

import sunpy
from sunpy.net import hek
from sunpy.net import attr


@pytest.fixture
def foostrwrap(request):
    return hek.attrs._StringParamAttrWrapper("foo")


@pytest.fixture
@pytest.mark.remote_data
def hek_client_creator():
    startTime = '2011/08/09 07:23:56'
    endTime = '2011/08/09 12:40:29'
    eventType = 'FL'

    hekTime = hek.attrs.Time(startTime, endTime)
    hekEvent = hek.attrs.EventType(eventType)

    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    return hek_query

def test_eventtype_collide():
    with pytest.raises(TypeError):
        hek.attrs.AR & hek.attrs.CE
    with pytest.raises(TypeError):
        (hek.attrs.AR & hek.attrs.Time((2011, 1, 1),
                                       (2011, 1, 2))) & hek.attrs.CE
        with pytest.raises(TypeError):
            (hek.attrs.AR | hek.attrs.Time((2011, 1, 1),
                                           (2011, 1, 2))) & hek.attrs.CE


def test_eventtype_or():
    assert (hek.attrs.AR | hek.attrs.CE).item == "ar,ce"


def test_paramattr():
    res = hek.attrs.walker.create(hek.attrs._ParamAttr("foo", "=", "bar"), {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '=', 'param0': 'foo'}


def test_stringwrapper_eq(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap == "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '=', 'param0': 'foo'}


def test_stringwrapper_lt(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap < "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '<', 'param0': 'foo'}


def test_stringwrapper_gt(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap > "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '>', 'param0': 'foo'}


def test_stringwrapper_le(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap <= "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '<=', 'param0': 'foo'}


def test_stringwrapper_ge(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap >= "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '>=', 'param0': 'foo'}


def test_stringwrapper_ne(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap != "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '!=', 'param0': 'foo'}


def test_stringwrapper_like(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap.like("bar"), {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': 'like', 'param0': 'foo'}


def test_err_dummyattr_create():
    with pytest.raises(TypeError):
        hek.attrs.walker.create(attr.DummyAttr(), {})


def test_err_dummyattr_apply():
    with pytest.raises(TypeError):
        hek.attrs.walker.apply(attr.DummyAttr(), {})


@pytest.mark.remote_data
def test_hek_client():
    startTime = '2011/08/09 07:23:56'
    endTime = '2011/08/09 12:40:29'
    eventType = 'FL'

    hekTime = hek.attrs.Time(startTime, endTime)
    hekEvent = hek.attrs.EventType(eventType)

    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    assert type(hek_query) == sunpy.net.hek.hek.HEKTable


@pytest.mark.remote_data
def test_hek_empty_search_result():
    startTime = '1985-05-04 00:00:00'
    endTime = '1985-05-04 00:00:00'
    eventType = 'FL'

    hekTime = hek.attrs.Time(startTime, endTime)
    hekEvent = hek.attrs.EventType(eventType)

    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    assert type(hek_query) == sunpy.net.hek.hek.HEKTable
    assert len(hek_query) == 0


@pytest.mark.remote_data
def test_getitem(hek_client_creator):
    hc = hek_client_creator
    assert hc.__getitem__(0) == hc[0]


@pytest.mark.remote_data
def test_get_voevent(hek_client_creator):
    hc = hek_client_creator
    ve = hc[0].get_voevent()
    assert len(ve['voe:VOEvent']) == 7


@pytest.mark.remote_data
def test_vso_time(hek_client_creator):
    hc = hek_client_creator
    ve = hc[0].vso_time
    assert type(ve) == sunpy.net.vso.attrs.Time


@pytest.mark.remote_data
def test_vso_instrument(hek_client_creator):
    hc = hek_client_creator
    try:
        vc = hc[1].vso_instrument
        assert type(vc) == sunpy.net.vso.attrs.Instrument
    except ValueError:
        assert 1


@pytest.mark.remote_data
def test_HEKRow_get(hek_client_creator):
    hc = hek_client_creator
    assert hc[0]['event_peaktime'] == hc[0].get('event_peaktime')
    assert hc[0].get('') is None
