# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from urllib2 import urlopen
from urllib import urlencode

from sunpy.net import attr


class ParamAttr(attr.ValueAttr):
    def __init__(self, name, op, value):
        attr.ValueAttr.__init__(self, [name])
        self.name = name
        self.op = op
        self.value = value


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


if __name__ == '__main__':
    attrs = ParamAttr('FRM_Name', '=', 'Karel Schrijver') & ParamAttr('FL_GOESCls', '>', 'X')
    print urlencode(walker.create(attrs, [0])[0])
