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


class AttrWalker(object):
    def __init__(self, appliers=None, creators=None, converters=None):
        self.methods = defaultdict(dict)
        self.cache = defaultdict(dict)
    
    def get_item(self, obj, method, super_=False):
        ocls = obj.__class__
        
        cached = self.cache[method].get(ocls, None)
        if cached is not None:
            return cached
        
        dct = self.methods[method]
        for cls in ocls.__mro__[super_:]:
            meth = dct.get(cls, None)
            if meth is not None:
                self.cache[method][ocls] = meth
                return meth
        raise KeyError
    
    def add(self, method, *types):
        self.cache[method] = {}
        def _dec(fun):
            for type_ in types:
                self.methods[method][type_] = fun
            return fun
        return _dec
    
    def add_creator(self, *types):
        return self.add('create', *types)
    
    def add_applier(self, *types):
        return self.add('apply', *types)
    
    def add_converter(self, *types):
        def _dec(fun):
            for type_ in types:
                self.converters[type_] = fun
            return fun
        return _dec
    
    def call(self, method, root, *args, **kwargs):
        return self.get_item(root, method)(self, root, *args, **kwargs)
    
    def super(self, method, root, *args, **kwargs):
        return self.get_item(root, method, True)(
            self, root, *args, **kwargs)
    
    def create(self, *args, **kwargs):
        return self.call('create', *args, **kwargs)

    def apply(self, *args, **kwargs):
        return self.call('apply', *args, **kwargs)
    
    def super_create(self, root, *args, **kwargs):
        return self.super('create', *args, **kwargs)

    def super_apply(self, root, *args, **kwargs):
        return self.super('apply', *args, **kwargs)
    
    def __copy__(self):
        return AttrWalker(self.appliers.copy(), self.creators.copy())


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
