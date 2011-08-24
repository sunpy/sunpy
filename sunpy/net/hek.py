# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

import sys

from urllib2 import urlopen
from urllib import urlencode
from datetime import datetime
from functools import partial

from sunpy.net import attr

DEFAULT_URL = 'http://www.lmsal.com/hek/her'


class ParamAttr(attr.ValueAttr):
    def __init__(self, name, op, value):
        attr.ValueAttr.__init__(self, {(name, op) : value})
        self.name = name
        self.op = op
        self.value = value


class BoolParamAttr(ParamAttr):
    def __init__(self, name, value='true'):
        ParamAttr.__init__(self, name, '=', value)
    
    def __neg__(self):
        if self.value == 'true':
            return BoolParamAttr(self.name, 'false')
        else:
            return BoolParamAttr(self.name)
    
    def __pos__(self):
        return BoolParamAttr(self.name)


class ListAttr(attr.Attr):
    def __init__(self, key, item):
        attr.Attr.__init__(self)
        
        self.key = key
        self.item = item
    
    def collides(self, other):
        return False
    
    def __eq__(self, other):
        if not isinstance(other, ListAttr):
            return False
        return self.key == other.key and self.item == other.item
    
    def __hash__(self):
        return hash((self.key, self.item))



class Time(attr.Attr):
    def __init__(self, start, end):
        attr.Attr.__init__(self)
        self.start = start
        self.end = end
    
    def collides(self, other):
        return isinstance(other, Time)
    
    def __eq__(self, other):
        if not isinstance(other, ListAttr):
            return False
        return self.start == other.start and self.end == other.end
    
    def __hash__(self):
        return hash((self.start, self.end))
    
    @classmethod
    def dt(cls, start, end):
        return cls(datetime(*start), datetime(*end))


class SpartialRegion(attr.Attr):
    def __init__(
        self, x1=-1200, y1=-1200, x2=1200, y2=1200, sys='helioprojective'):
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.sys = sys
    
    def collides(self, other):
        return isinstance(other, SpartialRegion)


walker = attr.AttrWalker()

@walker.add_creator(Time, SpartialRegion, ListAttr, ParamAttr, attr.AttrAnd)
def _c(walker, root, id_):
    value = {}
    walker.apply(root, id_, value)
    return [value]

@walker.add_applier(Time)
def _a(walker, root, id_, dct):
    dct['event_starttime'] = root.start.strftime('%Y-%m-%dT%H:%M:%S')
    dct['event_endtime'] = root.end.strftime('%Y-%m-%dT%H:%M:%S')
    return dct

@walker.add_applier(SpartialRegion)
def _a(walker, root, id_, dct):
    dct['x1'] = root.x1
    dct['y1'] = root.y1
    dct['x2'] = root.x2
    dct['y2'] = root.y2
    dct['event_coordsys'] = root.sys
    return dct

@walker.add_applier(ListAttr)
def _a(walker, root, id_, dct):
    if root.key in dct:
        dct[root.key] += ',%s' % root.item
    else:
        dct[root.key] = root.item
    return dct

@walker.add_applier(ParamAttr)
def _a(walker, root, id_, dct):
    nid = id_[0]
    dct['param%d' % nid] = root.name
    dct['op%d' % nid] = root.op
    dct['value%d' % nid] = root.value
    id_[0] += 1
    return dct

@walker.add_applier(attr.AttrAnd)
def _a(walker, root, id_, dct):
    for attr in root.attrs:
        walker.apply(attr, id_, dct)

@walker.add_creator(attr.AttrOr)
def _c(walker, root, id_):
    blocks = []
    for attr in self.attrs:
        blocks.extend(walker.create(attr, id_))
    return blocks

@walker.add_creator(attr.DummyAttr)
def _c(walker, root, id_):
    return {}

@walker.add_applier(attr.DummyAttr)
def _a(walker, root, id_, dct):
    pass


class ComparisonParamAttrWrapper(object):
    def __init__(self, name):
        self.name = name
    
    def __lt__(self, other):
        return ParamAttr(self.name, '<', other)
    
    def __le__(self, other):
        return ParamAttr(self.name, '<=', other)
    
    def __gt__(self, other):
        return ParamAttr(self.name, '>', other)
    
    def __ge__(self, other):
        return ParamAttr(self.name, '>=', other)
    
    def __eq__(self, other):
        return ParamAttr(self.name, '=', other)
    
    def __neq__(self, other):
        return ParamAttr(self.name, '!=', other)


class StringParamAttrWrapper(ComparisonParamAttrWrapper):
    def __init__(self, name):
        self.name = name
    
    def like(self, other):
        return ParamAttr(self.name, 'like', other)


class NumberParamAttrWrapper(ComparisonParamAttrWrapper):
    pass

EVENTS = [
    'AR', 'CE', 'CD', 'CH', 'CW', 'FI', 'FE', 'FA', 'FL', 'LP', 'OS', 'SS',
    'EF', 'CJ', 'PG', 'OT', 'NR', 'SG', 'SP', 'CR', 'CC', 'ER', 'TO'
]

    
class HEKClient(object):
    fields = {
        'FRM_Name': StringParamAttrWrapper,
        'FRM_HUMANFLAG': BoolParamAttr,
    }
    
    for elem in EVENTS:
        fields[elem] = partial(ListAttr, 'event_type')
    
    default = {
        'cosec': '2',
        'cmd': 'search',
        'type': 'column',
        'event_type': '**',
    }
    # Default to full disk.
    walker.apply(SpartialRegion(), [], default)
    
    def __init__(self, url=DEFAULT_URL):
        self.url = url
    
    def query(self, *query):
        if len(query) > 1:
            query = attr.and_(*query)
        
        data = walker.create(query, [0])
        ndata = []
        for elem in data:
            new = self.default.copy()
            new.update(elem)
            ndata.append(new)
        
        return urlopen(self.url, urlencode(ndata[0]))
    
    def __getitem__(self, name):
        return self.fields[name](name)


if __name__ == '__main__':
    import json
    import pprint
    c = HEKClient()
    pprint.pprint(json.load(c.query(
        Time.dt((2010, 1, 1), (2010, 1, 1, 1)),
        c['AR'],
    )))
