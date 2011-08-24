# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from urllib2 import urlopen
from urllib import urlencode

from sunpy.net import attr


class ParamAttr(attr.ValueAttr):
    def __init__(self, name, op, value):
        attr.ValueAttr.__init__(self, [(name, op)])
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


walker = attr.AttrWalker()

@walker.add_creator(ParamAttr)
def _c(walker, root, id_):
    value = {}
    walker.apply(root, id_, value)
    return [value]

@walker.add_applier(ParamAttr)
def _a(walker, root, id_, dct):
    nid = id_[0]
    dct['param%d' % nid] = root.name
    dct['op%d' % nid] = root.op
    dct['value%d' % nid] = root.value
    id_[0] += 1
    return dct

@walker.add_creator(attr.AttrAnd)
def _c(walker, root, id_):
    value = {}
    walker.apply(root, id_, value)
    return [value]

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


class StringParamAttrWrapper(object):
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
    
    def like(self, other):
        return ParamAttr(self.name, 'like', other)


class NumberParamAttrWrapper(object):
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


if __name__ == '__main__':
    FRM_Name = StringParamAttrWrapper('FRM_Name')
    FRM_HUMANFLAG = BoolParamAttr('FRM_HUMANFLAG')
    
    print walker.create(FRM_Name.like('Hello%'), [0])[0]
    print walker.create(FRM_Name < 'Hello', [0])[0]
    
    attrs = attr.and_(FRM_Name < 'Hello', -FRM_HUMANFLAG)
    print attrs
    print urlencode(walker.create(attrs, [0])[0])
