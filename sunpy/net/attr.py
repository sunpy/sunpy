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
        self.appliers = {} if appliers is None else appliers
        self.creators = {} if creators is None else creators
        self.converters = {} if converters is None else converters
    
    def get_item(self, obj, dct, super_=False):
        for cls in obj.__class__.__mro__[super_:]:
            if cls in self.converters:
                obj = dct[cls]()
                break
        
        for cls in obj.__class__.__mro__[super_:]:
            if cls in dct:
                return dct[cls]
        raise KeyError
    
    def add_creator(self, *types):
        def _dec(fun):
            for type_ in types:
                self.creators[type_] = fun
            return fun
        return _dec
    
    def add_applier(self, *types):
        def _dec(fun):
            for type_ in types:
                self.appliers[type_] = fun
            return fun
        return _dec
    
    def add_converter(self, *types):
        def _dec(fun):
            for type_ in types:
                self.converters[type_] = fun
            return fun
        return _dec
    
    def create(self, root, *args, **kwargs):
        return self.get_item(root, self.creators)(self, root, *args, **kwargs)

    def apply(self, root, *args, **kwargs):
        return self.get_item(root, self.appliers)(self, root, *args, **kwargs)
    
    def super_create(self, root, *args, **kwargs):
        return self.get_item(root, self.creators, True)(
            self, root, *args, **kwargs)

    def super_apply(self, root, *args, **kwargs):
        return self.get_item(root, self.appliers, True)(
            self, root, *args, **kwargs)
    
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
