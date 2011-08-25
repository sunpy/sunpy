from collections import defaultdict

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
        new = _DummyAttr()
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


class ValueAttr(Attr):
    def __init__(self, attrs):
        self.attrs = attrs
    
    def __repr__(self):
        return "<ValueAttr(%r)>" % (self.attrs)
    
    def __eq__(self, other):
        if not isinstance(other, ValueAttr):
            return False
        return self.attrs == other.attrs
    
    def __hash__(self):
        return hash(frozenset(self.attrs.iteritems()))
    
    def collides(self, other):
        if not isinstance(other, ValueAttr):
            return False
        return any(k in other.attrs for k in self.attrs)


class MultiMethod(object):
    def __init__(self, n=None, key=None):
        if n is not None and key is not None:
            raise ValueError
        elif n is None and key is None:
            raise ValueError
        
        self.n = n
        self.key = key
        
        self.methods = {}
        self.cache = {}
    
    def add(self, *types):
        self.cache = {}
        def _dec(fun):
            for type_ in types:
                self.methods[type_] = fun
            return fun
        return _dec
    
    def get_item(self, obj, super_=False):
        ocls = obj.__class__
        
        cached = self.cache.get(ocls, None)
        if cached is not None:
            return cached
        
        dct = self.methods
        for cls in ocls.__mro__[super_:]:
            meth = dct.get(cls, None)
            if meth is not None:
                self.cache[ocls] = meth
                return meth
        raise KeyError
    
    def __call__(self, *args, **kwargs):
        if self.n is not None:
            obj = args[self.n]
        else:
            obj = kwargs[self.key]
        return self.get_item(obj)(*args, **kwargs)
    
    def super(self, *args, **kwargs):
        if self.n is not None:
            obj = args[self.n]
        else:
            obj = kwargs[self.key]
        return self.get_item(obj, True)(*args, **kwargs)
        
    
class AttrWalker(object):
    def __init__(self):
        self.applymm = MultiMethod(1)
        self.createmm = MultiMethod(1)
    
    def add_creator(self, *types):
        return self.createmm.add(*types)
    
    def add_applier(self, *types):
        return self.applymm.add(*types)
    
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
