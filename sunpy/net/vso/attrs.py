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

from datetime import datetime

from astropy import units as u

from sunpy.time import TimeRange as _TimeRange
from sunpy.net.attr import (
    Attr, AttrWalker, AttrAnd, AttrOr, DummyAttr, ValueAttr
)
from sunpy.util.multimethod import MultiMethod
from sunpy.time import parse_time

__all__ = ['Wave', 'Time', 'Extent', 'Field', 'Provider', 'Source',
           'Instrument', 'Physobs', 'Pixels', 'Level', 'Resolution',
           'Detector', 'Filter', 'Sample', 'Quicklook', 'PScale']

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
    def __init__(self, wavemin, wavemax):
        if not all(isinstance(var, u.Quantity) for var in [wavemin, wavemax]):
            raise TypeError("Wave inputs must be astropy Quantities")

        # VSO just accept inputs as Angstroms, kHz or keV, the following
        # converts to any of these units depending on the spectral inputs
        # Note: the website asks for GHz, however it seems that using GHz
        # produces weird responses on VSO.
        convert = {'m': u.AA, 'Hz': u.kHz, 'eV': u.keV}
        for k in convert.keys():
            if wavemin.decompose().unit == (1 * u.Unit(k)).decompose().unit:
                unit = convert[k]
        try:
            self.min, self.max = sorted(
                value.to(unit) for value in [wavemin, wavemax]
                )
            self.unit = unit
        except NameError:
            raise ValueError("'{0}' is not a spectral supported unit".format(wavemin.unit))
        Attr.__init__(self)
        _Range.__init__(self, self.min, self.max, self.__class__)

    def collides(self, other):
        return isinstance(other, self.__class__)

    def __repr__(self):
        return "<Wave({0!r}, {1!r}, '{2!s}')>".format(self.min.value,
                                                      self.max.value,
                                                      self.unit)


class Time(Attr, _Range):
    """
    Specify the time range of the query.

    Parameters
    ----------

    start : SunPy Time String or `~sunpy.time.TimeRange`.
        The start time in a format parseable by `~sunpy.time.parse_time` or
        a `sunpy.time.TimeRange` object.

    end : SunPy Time String
        The end time of the range.

    near: SunPy Time String
        Return a singular record closest in time to this value as possible,
        inside the start and end window. Note: not all providers support this.

    """
    def __init__(self, start, end=None, near=None):
        if end is None and not isinstance(start, _TimeRange):
            raise ValueError("Specify start and end or start has to be a TimeRange")
        if isinstance(start, _TimeRange):
            self.start = start.start
            self.end = start.end
        else:
            self.start = parse_time(start)
            self.end = parse_time(end)
        self.near = None if near is None else parse_time(near)

        _Range.__init__(self, self.start, self.end, self.__class__)
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
        return '<Time({s.start!r}, {s.end!r}, {s.near!r})>'.format(s=self)

class Extent(Attr):
    # pylint: disable=R0913
    def __init__(self, x, y, width, length, atype):
        Attr.__init__(self)

        self.x = x
        self.y = y
        self.width = width
        self.length = length
        self.type = atype

    def collides(self, other):
        return isinstance(other, self.__class__)


class Field(ValueAttr):
    def __init__(self, fielditem):
        ValueAttr.__init__(self, {
            ('field', 'fielditem'): fielditem
        })


class _VSOSimpleAttr(Attr):
    """ A _SimpleAttr is an attribute that is not composite, i.e. that only
    has a single value, such as, e.g., Instrument('eit'). """
    def __init__(self, value):
        Attr.__init__(self)

        self.value = value

    def collides(self, other):
        return isinstance(other, self.__class__)

    def __repr__(self):
        return "<{cname!s}({val!r})>".format(
            cname=self.__class__.__name__, val=self.value)


class Provider(_VSOSimpleAttr):
    pass


class Source(_VSOSimpleAttr):
    pass


class Instrument(_VSOSimpleAttr):
    pass


class Physobs(_VSOSimpleAttr):
    pass


class Pixels(_VSOSimpleAttr):
    pass


class Level(_VSOSimpleAttr):
    pass


class Resolution(_VSOSimpleAttr):
    pass


class Detector(_VSOSimpleAttr):
    pass


class Filter(_VSOSimpleAttr):
    pass


class Sample(_VSOSimpleAttr):
    """
    Time interval for data sampling.

    Parameters
    ----------

    value : `astropy.units.Quantity`
        A sampling rate convertible to seconds.
    """
    @u.quantity_input(value=u.s)
    def __init__(self, value):
        super(Sample, self).__init__(value)
        self.value = value.to(u.s).value


class Quicklook(_VSOSimpleAttr):
    pass


class PScale(_VSOSimpleAttr):
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

walker.add_converter(_VSOSimpleAttr)(
    lambda x: ValueAttr({(x.__class__.__name__.lower(), ): x.value})
)

walker.add_converter(Wave)(
    lambda x: ValueAttr({
            ('wave', 'wavemin'): x.min.value,
            ('wave', 'wavemax'): x.max.value,
            ('wave', 'waveunit'): x.unit,
    })
)

# The idea of using a multi-method here - that means a method which dispatches
# by type but is not attached to said class - is that the attribute classes are
# designed to be used not only in the context of VSO but also elsewhere (which
# AttrAnd and AttrOr obviously are - in the HEK module). If we defined the
# filter method as a member of the attribute classes, we could only filter
# one type of data (that is, VSO data).
filter_results = MultiMethod(lambda *a, **kw: (a[0], ))

# If we filter with ANDed together attributes, the only items are the ones
# that match all of them - this is implementing  by ANDing the pool of items
# with the matched items - only the ones that match everything are there
# after this.
@filter_results.add_dec(AttrAnd)
def _(attr, results):
    res = set(results)
    for elem in attr.attrs:
        res &= filter_results(elem, res)
    return res

# If we filter with ORed attributes, the only attributes that should be
# removed are the ones that match none of them. That's why we build up the
# resulting set by ORing all the matching items.
@filter_results.add_dec(AttrOr)
def _(attr, results):
    res = set()
    for elem in attr.attrs:
        res |= filter_results(elem, results)
    return res

# Filter out items by comparing attributes.
@filter_results.add_dec(_VSOSimpleAttr)
def _(attr, results):
    attrname = attr.__class__.__name__.lower()
    return set(
        item for item in results
        # Some servers seem to omit some fields. No way to filter there.
        if not hasattr(item, attrname) or
        getattr(item, attrname).lower() == attr.value.lower()
    )

# The dummy attribute does not filter at all.
@filter_results.add_dec(DummyAttr, Field)
def _(attr, results):
    return set(results)


@filter_results.add_dec(Wave)
def _(attr, results):
    return set(
        it for it in results
        if
        it.wave.wavemax is not None
        and
        attr.min <= it.wave.wavemax.to(u.angstrom, equivalencies=u.spectral())
        and
        it.wave.wavemin is not None
        and
        attr.max >= it.wave.wavemin.to(u.angstrom, equivalencies=u.spectral())
    )

@filter_results.add_dec(Time)
def _(attr, results):
    return set(
        it for it in results
        if
        it.time.end is not None
        and
        attr.min <= datetime.strptime(it.time.end, TIMEFORMAT)
        and
        it.time.start is not None
        and
        attr.max >= datetime.strptime(it.time.start, TIMEFORMAT)
    )

@filter_results.add_dec(Extent)
def _(attr, results):
    return set(
        it for it in results
        if
        it.extent.type is not None
        and
        it.extent.type.lower() == attr.type.lower()
    )
