# pylint: disable=C0103,R0903

from datetime import datetime

from sunpy.net.attr import (
    Attr, ValueAttr, AttrWalker, AttrAnd, AttrOr, DummyAttr, and_, ValueAttr
)
from sunpy.util.util import to_angstrom

TIMEFORMAT = '%Y%m%d%H%M%S'

class Range(object):
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


class Wave(Attr, Range): 
    def __init__(self, wavemin, wavemax, waveunit='Angstrom'):        
        self.min, self.max = sorted(
            to_angstrom(v, waveunit) for v in [wavemin, wavemax]
        )
        self.unit = 'Angstrom'
        
        Attr.__init__(self)
        Range.__init__(self, self.min, self.max, self.__class__)
    
    def collides(self, other):
        return isinstance(other, self.__class__)


class Time(Attr, Range):
    def __init__(self, start, end, near=None):
        self.start = start
        self.end = end
        self.near = near

        Range.__init__(self, start, end, self.__class__)
        Attr.__init__(self)
    
    def collides(self, other):
        return isinstance(other, self.__class__)
    
    def __xor__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError
        if self.near is not None or other.near is not None:
            raise TypeError
        return Range.__xor__(self, other)
    
    @classmethod
    def dt(cls, start, end, near=None):
        if near is not None:
            near = datetime(*near)
        return cls(datetime(*start), datetime(*end), near)
    
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


class SimpleAttr(Attr):
    def __init__(self, value):
        Attr.__init__(self)
        
        self.value = value
    
    def collides(self, other):
        return isinstance(other, self.__class__)


class Provider(SimpleAttr):
    pass


class Source(SimpleAttr):
    pass


class Instrument(SimpleAttr):
    pass


class Physobj(SimpleAttr):
    pass


class Pixels(SimpleAttr):
    pass


class Level(SimpleAttr):
    pass


class Resolution(SimpleAttr):
    pass


class Detector(SimpleAttr):
    pass


class Filter(SimpleAttr):
    pass


class Sample(SimpleAttr):
    pass


class Quicklook(SimpleAttr):
    pass


class PScale(SimpleAttr):
    pass


# The walker specifies how the Attr-tree is converted to a query the
# server can handle.
walker = AttrWalker()

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
                x.near.strftime(TIMEFORMAT) if x.near is not None else ''),
    })
)

walker.add_converter(SimpleAttr)(
    lambda x: ValueAttr({(x.__class__.__name__.lower(), ): x.value})
)
