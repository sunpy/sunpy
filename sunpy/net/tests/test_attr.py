# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

from sunpy.net import attr
from sunpy.net.vso import attrs

def test_dummyattr():
    one = attr.DummyAttr()
    other = attr.ValueAttr({'a': 'b'})
    assert (one | other) is other
    assert (one & other) is other

def test_and_nesting():
    a = attr.and_(attrs.Level(0),
                  attr.AttrAnd((attrs.Instrument('EVE'),
                                attrs.Time("2012/1/1", "2012/01/02"))))
    # Test that the nesting has been removed.
    assert len(a.attrs) == 3

def test_or_nesting():
    a = attr.or_(attrs.Instrument('a'),
                  attr.AttrOr((attrs.Instrument('b'),
                               attrs.Instrument('c'))))
    # Test that the nesting has been removed.
    assert len(a.attrs) == 3
