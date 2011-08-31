from itertools import chain, repeat

from sunpy.util.multimethod import MultiMethod

class Attr(object):
    def __and__(self, other):
        if isinstance(other, AttrOr):
            return AttrOr([elem & self for elem in other.attrs])
        if self.collides(other):
            return NotImplemented
        return AttrAnd([self, other])
    
    def __or__(self, other):
        # Optimization.
        if self == other:
            return self
        return AttrOr([self, other])
    
    def collides(self, other):
        raise NotImplementedError


class DummyAttr(Attr):
    def __and__(self, other):
        return other
    
    def __or__(self, other):
        return other
    
    def collides(self, other):
        return False
    
    def __hash__(self):
        return hash(None)
    
    def __eq__(self, other):
        return isinstance(other, DummyAttr)


class AttrAnd(Attr):
    def __init__(self, attrs):
        self.attrs = attrs
    
    def __and__(self, other):
        if any(other.collides(elem) for elem in self.attrs):
            return NotImplemented
        if isinstance(other, AttrAnd):
            return AttrAnd(self.attrs + other.attrs)
        if isinstance(other, AttrOr):
            return AttrOr([elem & self for elem in other.attrs])
        return AttrAnd(self.attrs + [other])
    
    __rand__ = __and__
    
    def __repr__(self):
        return "<AttrAnd(%r)>" % self.attrs
    
    def __eq__(self, other):
        if not isinstance(other, AttrAnd):
            return False
        return set(self.attrs) == set(other.attrs)
    
    def __hash__(self):
        return hash(frozenset(self.attrs))
    
    def collides(self, other):
        return any(elem.collides(other) for elem in self)


class AttrOr(Attr):
    def __init__(self, attrs):
        self.attrs = attrs
    
    def __or__(self, other):
        if isinstance(other, AttrOr):
            return AttrOr(self.attrs + other.attrs)
        return AttrOr(self.attrs + [other])
    
    __ror__ = __or__
    
    def __and__(self, other):
        return AttrOr([elem & other for elem in self.attrs])
    
    __rand__ = __and__
    
    def __xor__(self, other):
        new = AttrOr([])
        for elem in self.attrs:
            new |= elem ^ other
        return new
    
    def __repr__(self):
        return "<AttrOr(%r)>" % self.attrs
    
    def __eq__(self, other):
        if not isinstance(other, AttrOr):
            return False
        return set(self.attrs) == set(other.attrs)
    
    def __hash__(self):
        return hash(frozenset(self.attrs))
    
    def collides(self, other):
        return all(elem.collides(other) for elem in self)


class KeysAttr(Attr):
    def __init__(self, attrs):
        self.attrs = attrs
    
    def __repr__(self):
        return "<KeysAttr(%r)>" % (self.attrs)
    
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.attrs == other.attrs
    
    def __hash__(self):
        return hash(frozenset(self.attrs))
    
    def collides(self, other):
        if not isinstance(other, ValueAttr):
            return False
        return any(k in other.attrs for k in self.attrs)


class ValueAttr(KeysAttr):
    def __init__(self, attrs):
        KeysAttr.__init__(self, attrs)
        
        self.attrs = attrs
    
    def __repr__(self):
        return "<ValueAttr(%r)>" % (self.attrs)

    def __hash__(self):
        return hash(frozenset(self.attrs.iteritems()))

    
class AttrWalker(object):
    def __init__(self):
        self.applymm = MultiMethod(lambda *a, **kw: (a[1], ))
        self.createmm = MultiMethod(lambda *a, **kw: (a[1], ))
    
    def add_creator(self, *types):
        def _dec(fun):
            for type_ in types:
                self.createmm.add(fun, (type_, ))
            return fun
        return _dec
    
    def add_applier(self, *types):
        def _dec(fun):
            for type_ in types:
                self.applymm.add(fun, (type_, ))
            return fun
        return _dec
    
    def add_converter(self, *types):
        def _dec(fun):
            for type_ in types:                
                self.applymm.add(self.cv_apply(fun), (type_, ))
                self.createmm.add(self.cv_create(fun), (type_, ))
            return fun
        return _dec
    
    def cv_apply(self, fun):
        def _fun(*args, **kwargs):
            nargs, nkwargs = fun(*args, **kwargs)
            return self.apply(*nargs, **nkwargs)
        return _fun
    
    def cv_create(self, fun):
        def _fun(*args, **kwargs):
            nargs, nkwargs = fun(*args, **kwargs)
            return self.create(*nargs, **nkwargs)
        return _fun
    
    def create(self, *args, **kwargs):
        return self.createmm(self, *args, **kwargs)

    def apply(self, *args, **kwargs):
        return self.applymm(self, *args, **kwargs)
    
    def super_create(self, *args, **kwargs):
        return self.createmm.super(self, *args, **kwargs)

    def super_apply(self, *args, **kwargs):
        return self.applymm.super(self, *args, **kwargs)


def and_(*args):
    """ Trick operator precendence. 
    
    and_(foo < bar, bar < baz)
    """
    value = DummyAttr()
    for elem in args:
        value &= elem
    return value

def or_(*args):
    """ Trick operator precendence. 
    
    or_(foo < bar, bar < baz)
    """
    value = DummyAttr()
    for elem in args:
        value |= elem
    return value
