# -*- coding: utf-8 -*-
"""
Allow representation of queries as logic expressions. This module makes
sure that attributes that are combined using the two logic operations AND (&)
and OR (|) always are in disjunctive normal form, that is, there are only two
levels ­- the first being disjunction and the second being conjunction. In other
words, every combinations of attributes looks like this:
(a AND b AND c) OR (d AND e).

Walkers are used to traverse the tree that results from combining attributes.
They are implemented using sunpy.util.multimethod. Multimethods are functions
that are not assigned to classes but still dispatch by type of one or more
of their arguments. For more information about multimethods, refer to
sunpy.util.multimethod.

Please note that & is evaluated first, so A & B | C is equivalent to
(A & B) | C.
"""
import keyword
import warnings
from collections import defaultdict, namedtuple

from astropy.utils.misc import isiterable

from sunpy.util.multimethod import MultiMethod
from sunpy.extern.six import iteritems

_ATTR_TUPLE = namedtuple("attr", "name name_long desc")

__all__ = ['Attr', 'DummyAttr', 'SimpleAttr', 'AttrAnd',
           'AttrOr', 'ValueAttr', 'and_', 'or_']


def make_tuple():
    return _ATTR_TUPLE([], [], [])


class AttrMeta(type):
    """
    We want to enable automatic discovery via tab completion of subclasses of Attrs.
    To do this we have to create a metaclass that redefines the methods that Python uses to normally do this.
    This would allow for example that `attrs.Instrument` to be able to tab complete to `attrs.Instrument.aia`.
    However, it must be done by a downloader client which know what Attrs they support.
    """

    # The aim is to register Attrs as a namedtuple of lists
    # So we define the namedtuple outside of AttrMeta, see above.
    _attr_registry = defaultdict(make_tuple)

    def __getattr__(self, item):
        """
        Our method for Attrs is to register using `type(Attr)` as keys into a dictionary.
        _attr_registry is a dictionary with a key being `type(Attr)` and the value being the namedtuple of lists.
        As a result we index `_attr_registry` with `[self]` which will be the `type(Attr)` to access the dictionary.
        This will return the namedtuple that has three attributes: `name`, `name_long` and `desc`. Each of which are a list.
        `name` will be the attribute name, `name_long` is the original name passed in and `desc` the description of the object.
        """
        registry = self._attr_registry[self]  # Get the revelant entries.
        names = registry.name  # All the attribute names under that type(Attr)
        if item in names:
            return self(registry.name_long[names.index(item)])  # We return Attr(name_long) to create the Attr requested.
        else:
            raise AttributeError

    def __dir__(self):
        """
        To tab complete in Python we need to add to the `__dir__()` return.
        So we add the attribute name of any Attrs which have been added to the _attr_registry here.
        """
        return super().__dir__() + self._attr_registry[self].name

#    def __repr__(self):
#       """
#       This enables the pretty printing of Attrs.
#       """
#        return str()


class Attr(metaclass=AttrMeta):
    """ This is the base for all attributes. """
    def __and__(self, other):
        if isinstance(other, AttrOr):
            return AttrOr([elem & self for elem in other.attrs])
        if self.collides(other):
            return NotImplemented
        if isinstance(other, AttrAnd):
            return AttrAnd([self] + list(other.attrs))
        return AttrAnd([self, other])

    def __hash__(self):
        return hash(frozenset(iteritems(vars(self))))

    def __or__(self, other):
        # Optimization.
        if self == other:
            return self
        if isinstance(other, AttrOr):
            return AttrOr([self] + list(other.attrs))
        return AttrOr([self, other])

    def collides(self, other):
        raise NotImplementedError

    def __eq__(self, other):
        return dict(vars(self)) == dict(vars(other))

    @classmethod
    def update_values(cls, adict):
        """
        This is the method that clients will use to register their `Attrs`.
        The input should be a dictionary.
        Each key should be type(attrs.x) of any attrs require to exist for `_can_handle_query` to return successfully.
        Please note that you can pass in type(attrs.x(Input)) or just the attrs.x object.
        The value for each key should be a list of tuples.
        The tuple should be of the form (`Name`, `Description`).
        If you do not want to add a description, you can put None or "".
        We sanitize the name you provide.
        We remove all special characters and make it all lower case.
        If it still isn't valid we will append to the start of the name to make it a valid attribute name.

        We have an example here for an Instrument Class.

        Example
        -------
        >>> from sunpy.net import attr, attrs
        >>> attr.Attr.update_values({attrs.Instrument: [('AIA', 'AIA is in Space.'), ('HMI', 'HMI is next to AIA.')]})
        >>> attr.Attr._attr_registry[attrs.Instrument]
        attr(name=['aia', 'hmi'], name_long=['AIA', 'HMI'], desc=['AIA is in Space.', 'HMI is next to AIA.'])
        >>> attr.Attr._attr_registry[attrs.Instrument].name
        ['aia', 'hmi']
        >>> attr.Attr._attr_registry[attrs.Instrument].name_long
        ['AIA', 'HMI']
        >>> attr.Attr._attr_registry[attrs.Instrument].desc
        ['AIA is in Space.', 'HMI is next to AIA.']
        >>> attrs.Instrument.aia, attrs.Instrument.hmi
        (<Instrument('AIA')>, <Instrument('HMI')>)
        """
        for k, v in adict.items():
            if isiterable(v) and not isinstance(v, str):
                for pair in v:
                    if len(pair) != 2:
                        raise ValueError('Please check your input dictionary, \
                                          it appears the value for a key value is not length 2 (`Name`, `Description`)')
                    else:
                        # Sanitize name, we remove all special characters and make it all lower case
                        name = ''.join(char for char in pair[0] if char.isalnum()).lower()
                        if keyword.iskeyword(name) or not name.isidentifier():
                                name = 'A' + name
                                warnings.warn("Attribute name has been appended to make it a valid identifier.")
                        cls._attr_registry[k][0].append(name)
                        cls._attr_registry[k][1].append(pair[0])
                        cls._attr_registry[k][2].append(pair[1])
            else:
                raise ValueError('Please check your input dictionary, it appears the value for a key is not iterable or a string')


class DummyAttr(Attr):
    """ Empty attribute. Useful for building up queries. Returns other
    attribute when ORed or ANDed. It can be considered an empty query
    that you can use as an initial value if you want to build up your
    query in a loop.

    So, if we wanted an attr matching all the time intervals between the times
    stored as (from, to) tuples in a list, we could do.

    attr = DummyAttr()
    for from\_, to in times:
        attr |= Time(from\_, to)
    """
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


class SimpleAttr(Attr):
    """
    An attribute that only has a single value.

    This type of attribute is not a composite and has a single value such as
    ``Instrument('EIT')``.

    Parameters
    ----------
    value : `object`
       The value for the attribute to hold.
    """
    def __init__(self, value):
        Attr.__init__(self)
        self.value = value

    def collides(self, other):
        return isinstance(other, self.__class__)

    def __repr__(self):
        return "<{cname!s}({val!r})>".format(
            cname=self.__class__.__name__, val=self.value)


class AttrAnd(Attr):
    """ Attribute representing attributes ANDed together. """
    def __init__(self, attrs):
        Attr.__init__(self)
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
        return "<AttrAnd({att!r})>".format(att=self.attrs)

    def __eq__(self, other):
        if not isinstance(other, AttrAnd):
            return False
        return set(self.attrs) == set(other.attrs)

    def __hash__(self):
        return hash(frozenset(self.attrs))

    def collides(self, other):
        return any(elem.collides(other) for elem in self.attrs)


class AttrOr(Attr):
    """ Attribute representing attributes ORed together. """
    def __init__(self, attrs):
        Attr.__init__(self)
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
            try:
                new |= elem ^ other
            except TypeError:
                pass
        return new

    def __contains__(self, other):
        for elem in self.attrs:
            try:
                if other in elem:
                    return True
            except TypeError:
                pass
        return False

    def __repr__(self):
        return "<AttrOr({att!r})>".format(att=self.attrs)

    def __eq__(self, other):
        if not isinstance(other, AttrOr):
            return False
        return set(self.attrs) == set(other.attrs)

    def __hash__(self):
        return hash(frozenset(self.attrs))

    def collides(self, other):
        return all(elem.collides(other) for elem in self.attrs)


class ValueAttr(Attr):
    def __init__(self, attrs):
        Attr.__init__(self)
        self.attrs = attrs

    def __repr__(self):
        return "<ValueAttr({att!r})>".format(att=self.attrs)

    def __hash__(self):
        return hash(frozenset(self.attrs.iteritems.items()))

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.attrs == other.attrs

    def collides(self, other):
        if not isinstance(other, self.__class__):
            return False
        return any(k in other.attrs for k in self.attrs)


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
            args = list(args)
            args[1] = fun(args[1])
            return self.applymm(*args, **kwargs)
        return _fun

    def cv_create(self, fun):
        def _fun(*args, **kwargs):
            args = list(args)
            args[1] = fun(args[1])
            return self.createmm(*args, **kwargs)
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
    """ Trick operator precedence.

    and_(foo < bar, bar < baz)
    """
    value = DummyAttr()
    for elem in args:
        value &= elem
    return value


def or_(*args):
    """ Trick operator precedence.

    or_(foo < bar, bar < baz)
    """
    value = DummyAttr()
    for elem in args:
        value |= elem
    return value
