class Attr(object):
    def __and__(self, other):
        if isinstance(other, AttrOr):
            return AttrOr([elem & self for elem in other.attrs])
        if isinstance(other, self.__class__):
            # A record cannot match two different values
            # for the same attribute.
            # TODO: Error?
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
    
    def apply(self, queryblock):
        pass
    
    def create(self, api):
        return api.factory.create('QueryRequestBlock')
    
    def collides(self, other):
        return False


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
    
    def apply(self, queryblock):
        for attr in self.attrs:
            attr.apply(queryblock)
    
    def create(self, api):
        # TODO: Prove that we can assume that only _SimpleAttr and
        # _ComplexAttr can possibly exist here.
        value = api.factory.create('QueryRequestBlock')
        self.apply(value)
        return [value]
    
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
    
    def apply(self, queryblock):
        # TODO: Prove this is unreachable.
        raise NotImplementedError
    
    def create(self, api):
        blocks = []
        for attr in self.attrs:
            blocks.extend(attr.create(api))
        return blocks
    
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
    
    def apply(self, queryblock):       
        for k, v in self.attrs.iteritems():
            lst = k[-1]
            rest = k[:-1]
            
            block = queryblock
            for elem in rest:
                block = block[elem]
            block[lst] = v
    
    def create(self, api):
        value = api.factory.create('QueryRequestBlock')
        self.apply(value)
        return [value]
    
    def __repr__(self):
        return "<ValueAttr(%r, %r)>" % (self.attr, self.attrs)
    
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
