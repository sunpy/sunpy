# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>
#
# This module was developed with funding provided by
# the ESA Summer of Code (2011).
#
# pylint: disable=C0103,R0903

"""
Attributes that can be used to construct VSO queries. Attributes are the
fundamental building blocks of queries that, together with the two
operations of AND and OR (and in some rare cases XOR) can be used to
construct complex queries. Most attributes can only be used once in an
AND-expression, if you still attempt to do so it is called a collision,
for a quick example think about how the system should handle
Instrument('aia') & Instrument('eit').
"""

from __future__ import absolute_import

from sunpy.net.attr import (
    Attr, ValueAttr, AttrWalker, AttrAnd, AttrOr, DummyAttr, ValueAttr
)
from sunpy.util.util import to_angstrom
from sunpy.time import parse_time

TIMEFORMAT = '%Y%m%d%H%M%S'

class _Range(object):
    def __init__(self, min_, max_, create):
        self.min = min_
        self.max = max_
        self.create = create
    
    def __xor__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        
        new = DummyAttr()
        if self.min < other.min:            
            new |= self.create(self.min, min(other.min, self.max))
        if other.max < self.max:
            new |= self.create(other.max, self.max)
        return new
    
    def __contains__(self, other):
        return self.min <= other.min and self.max >= other.max


class Wave(Attr, _Range):
    def __init__(self, wavemin, wavemax, waveunit='Angstrom'):        
        self.min, self.max = sorted(
            to_angstrom(v, waveunit) for v in [wavemin, wavemax]
        )
        self.unit = 'Angstrom'
        
        Attr.__init__(self)
        _Range.__init__(self, self.min, self.max, self.__class__)
    
    def collides(self, other):
        return isinstance(other, self.__class__)


class Time(Attr, _Range):
    def __init__(self, start, end, near=None):
        self.start = parse_time(start)
        self.end = parse_time(end)
        self.near = None if near is None else parse_time(near)

        _Range.__init__(self, start, end, self.__class__)
        Attr.__init__(self)
    
    def collides(self, other):
        return isinstance(other, self.__class__)
    
    def __xor__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError
        if self.near is not None or other.near is not None:
            raise TypeError
        return _Range.__xor__(self, other)
    
    def pad(self, timedelta):
        return Time(self.start - timedelta, self.start + timedelta)
    
    def __repr__(self):
        return '<Time(%r, %r, %r)>' % (self.start, self.end, self.near)


class Extent(Attr):
    # pylint: disable=R0913
    def __init__(self, x, y, width, length, type_):
        Attr.__init__(self)
        
        self.x = x
        self.y = y
        self.width = width
        self.length = length
        self.type = type_
    
    def collides(self, other):
        return isinstance(other, self.__class__)


class Field(ValueAttr):
    def __init__(self, fielditem):
        ValueAttr.__init__(self, {
            ['field', 'fielditem']: fielditem
        })


class _SimpleAttr(Attr):
    def __init__(self, value):
        Attr.__init__(self)
        
        self.value = value
    
    def collides(self, other):
        return isinstance(other, self.__class__)
    
    def __repr__(self):
        return "<%s(%r)>" % (self.__class__.__name__, self.value)


class Provider(_SimpleAttr):
    pass


class Source(_SimpleAttr):
    pass


class Instrument(_SimpleAttr):
    pass


class Physobs(_SimpleAttr):
    pass


class Pixels(_SimpleAttr):
    pass


class Level(_SimpleAttr):
    pass


class Resolution(_SimpleAttr):
    pass


class Detector(_SimpleAttr):
    pass


class Filter(_SimpleAttr):
    pass


class Sample(_SimpleAttr):
    pass


class Quicklook(_SimpleAttr):
    pass


class PScale(_SimpleAttr):
    pass


# The walker specifies how the Attr-tree is converted to a query the
# server can handle.
walker = AttrWalker()

# The _create functions make a new VSO query from the attribute tree,
# the _apply functions take an existing query-block and update it according
# to the attribute tree passed in as root. Different attributes require
# different functions for conversion into query blocks.

@walker.add_creator(ValueAttr, AttrAnd)
# pylint: disable=E0102,C0103,W0613
def _create(wlk, root, api):
    """ Implementation detail. """
    value = api.factory.create('QueryRequestBlock')
    wlk.apply(root, api, value)
    return [value]

@walker.add_applier(ValueAttr)
# pylint: disable=E0102,C0103,W0613
def _apply(wlk, root, api, queryblock):
    """ Implementation detail. """
    for k, v in root.attrs.iteritems():
        lst = k[-1]
        rest = k[:-1]
        
        block = queryblock
        for elem in rest:
            block = block[elem]
        block[lst] = v

@walker.add_applier(AttrAnd)
# pylint: disable=E0102,C0103,W0613
def _apply(wlk, root, api, queryblock):
    """ Implementation detail. """
    for attr in root.attrs:
        wlk.apply(attr, api, queryblock)

@walker.add_creator(AttrOr)
# pylint: disable=E0102,C0103,W0613
def _create(wlk, root, api):
    """ Implementation detail. """
    blocks = []
    for attr in root.attrs:
        blocks.extend(wlk.create(attr, api))
    return blocks

@walker.add_creator(DummyAttr)
# pylint: disable=E0102,C0103,W0613
def _create(wlk, root, api):
    """ Implementation detail. """
    return api.factory.create('QueryRequestBlock')

@walker.add_applier(DummyAttr)
# pylint: disable=E0102,C0103,W0613
def _apply(wlk, root, api, queryblock):
    """ Implementation detail. """
    pass


# Converters take a type unknown to the walker and convert it into one
# known to it. All of those convert types into ValueAttrs, which are
# handled above by just assigning according to the keys and values of the
# attrs member.
walker.add_converter(Extent)(
    lambda x: ValueAttr(
        dict((('extent', k), v) for k, v in vars(x).iteritems())
    )
)

walker.add_converter(Time)(
    lambda x: ValueAttr({
            ('time', 'start'): x.start.strftime(TIMEFORMAT),
            ('time', 'end'): x.end.strftime(TIMEFORMAT) ,
            ('time', 'near'): (
                x.near.strftime(TIMEFORMAT) if x.near is not None else None),
    })
)

walker.add_converter(_SimpleAttr)(
    lambda x: ValueAttr({(x.__class__.__name__.lower(), ): x.value})
)

walker.add_converter(Wave)(
    lambda x: ValueAttr({
            ('wave', 'wavemin'): x.min,
            ('wave', 'wavemax'): x.max,
            ('wave', 'waveunit'): x.unit,
    })
)
