# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

#pylint: disable=W0613

from __future__ import absolute_import

import pytest

from sunpy.net import hek
from sunpy.net import attr


@pytest.fixture
def foostrwrap(request):
    return hek.attrs._StringParamAttrWrapper("foo")


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
