# The template can be found in tools/hektemplate.py

from datetime import datetime
from sunpy.net import attr
from sunpy.util.util import anytim

class _ParamAttr(attr.Attr):
    def __init__(self, name, op, value):
        attr.Attr.__init__(self)
        self.name = name
        self.op = op
        self.value = value
    
    def collides(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.op == other.op and self.name == other.name


class _BoolParamAttr(_ParamAttr):
    def __init__(self, name, value='true'):
        _ParamAttr.__init__(self, name, '=', value)
    
    def __neg__(self):
        if self.value == 'true':
            return _BoolParamAttr(self.name, 'false')
        else:
            return _BoolParamAttr(self.name)
    
    def __pos__(self):
        return _BoolParamAttr(self.name)


class _ListAttr(attr.Attr):
    def __init__(self, key, item):
        attr.Attr.__init__(self)
        
        self.key = key
        self.item = item
    
    def collides(self, other):
        return False
    
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return vars(self) == vars(other)
    
    def __hash__(self):
        return hash(tuple(vars(self).itervalues()))


class EventType(_ListAttr):
    def __init__(self, item):
        _ListAttr.__init__(self, 'event_type', item)


class Time(attr.Attr):
    def __init__(self, start, end):
        attr.Attr.__init__(self)
        self.start = start
        self.end = end
    
    def collides(self, other):
        return isinstance(other, Time)
    
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return vars(self) == vars(other)
    
    def __hash__(self):
        return hash(tuple(vars(self).itervalues()))
    
    @classmethod
    def dt(cls, start, end):
        return cls(datetime(*start), datetime(*end))


# pylint: disable=R0913
class SpartialRegion(attr.Attr):
    def __init__(
        self, x1=-1200, y1=-1200, x2=1200, y2=1200, sys='helioprojective'):
        attr.Attr.__init__(self)
        
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.sys = sys
    
    def collides(self, other):
        return isinstance(other, SpartialRegion)
    
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return vars(self) == vars(other)
    
    def __hash__(self):
        return hash(tuple(vars(self).itervalues()))


class Contains(attr.Attr):
    def __init__(self, *types):
        attr.Attr.__init__(self)
        self.types = types
    
    def collides(self, other):
        return False
    
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return vars(self) == vars(other)
    
    def __hash__(self):
        return hash(tuple(vars(self).itervalues()))


class ComparisonParamAttrWrapper(object):
    def __init__(self, name):
        self.name = name
    
    def __lt__(self, other):
        return _ParamAttr(self.name, '<', other)
    
    def __le__(self, other):
        return _ParamAttr(self.name, '<=', other)
    
    def __gt__(self, other):
        return _ParamAttr(self.name, '>', other)
    
    def __ge__(self, other):
        return _ParamAttr(self.name, '>=', other)
    
    def __eq__(self, other):
        return _ParamAttr(self.name, '=', other)
    
    def __neq__(self, other):
        return _ParamAttr(self.name, '!=', other)


class _StringParamAttrWrapper(ComparisonParamAttrWrapper):
    def like(self, other):
        return _ParamAttr(self.name, 'like', other)


class _NumberParamAttrWrapper(ComparisonParamAttrWrapper):
    pass


walker = attr.AttrWalker()

@walker.add_applier(Contains)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    dct['type'] = 'contains'
    if not Contains in state:
        state[Contains] = 1
    
    nid = state[Contains]
    n = 0
    for n, type_ in enumerate(root.types):
        dct['event_type%d' % (nid + n)] = type_
    state[Contains] += n
    return dct

@walker.add_creator(
    Time, SpartialRegion, _ListAttr, _ParamAttr, attr.AttrAnd, Contains)
# pylint: disable=E0102,C0103,W0613
def _c(wlk, root, state):
    value = {}
    wlk.apply(root, state, value)
    return [value]

@walker.add_applier(Time)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    dct['event_starttime'] = anytim(root.start).strftime('%Y-%m-%dT%H:%M:%S')
    dct['event_endtime'] = anytim(root.end).strftime('%Y-%m-%dT%H:%M:%S')
    return dct

@walker.add_applier(SpartialRegion)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    dct['x1'] = root.x1
    dct['y1'] = root.y1
    dct['x2'] = root.x2
    dct['y2'] = root.y2
    dct['event_coordsys'] = root.sys
    return dct

@walker.add_applier(_ListAttr)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    if root.key in dct:
        dct[root.key] += ',%s' % root.item
    else:
        dct[root.key] = root.item
    return dct

@walker.add_applier(EventType)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    if dct.get('type', None) == 'contains':
        raise ValueError
    
    return wlk.super_apply(super(EventType, root), state, dct)

@walker.add_applier(_ParamAttr)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    if not _ParamAttr in state:
        state[_ParamAttr] = 0
    
    nid = state[_ParamAttr]
    dct['param%d' % nid] = root.name
    dct['op%d' % nid] = root.op
    dct['value%d' % nid] = root.value
    state[_ParamAttr] += 1
    return dct

@walker.add_applier(attr.AttrAnd)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    for attribute in root.attrs:
        wlk.apply(attribute, state, dct)

@walker.add_creator(attr.AttrOr)
# pylint: disable=E0102,C0103,W0613
def _c(wlk, root, state):
    blocks = []
    for attribute in root.attrs:
        blocks.extend(wlk.create(attribute, state))
    return blocks

@walker.add_creator(attr.DummyAttr)
# pylint: disable=E0102,C0103,W0613
def _c(wlk, root, state):
    return {}

@walker.add_applier(attr.DummyAttr)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    pass
