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
import re
import keyword
import warnings
from collections import defaultdict, namedtuple

from astropy.utils.misc import isiterable

from sunpy.util.multimethod import MultiMethod

_ATTR_TUPLE = namedtuple("attr", "name name_long desc")
_REGEX = re.compile(r"^[\d]([^\d].*)?$")

__all__ = ['Attr', 'DummyAttr', 'SimpleAttr', 'AttrAnd',
           'AttrOr', 'ValueAttr', 'and_', 'or_']


def make_tuple():
    return _ATTR_TUPLE([], [], [])


def strtonum(value):
    """
    For numbers 0 to 9, return the number spelled out. Otherwise, return the
    number. This follows Associated Press style.  This always returns a string
    unless the value was not int-able, unlike the Django filter.
    Taken from: https://github.com/jmoiron/humanize/blob/master/humanize/number.py#L81
    """
    try:
        value = int(value)
    except (TypeError, ValueError):
        return value
    if not 0 <= value < 10:
        return str(value)
    return ('zero', 'one', 'two', 'three', 'four', 'five', 'six',
            'seven', 'eight', 'nine')[value]


class AttrMeta(type):
    """
    We want to enable discovery, by tab completion, of values for all subclasses of Attr.
    So have to create a metaclass that overloads the methods that Python uses, so that they work on the classes.
    This would allow that `attrs.Instrument` to be able to tab complete to `attrs.Instrument.aia`.
    """

    # The aim is to register Attrs as a namedtuple of lists
    # So we define the namedtuple outside of AttrMeta, see above.
    _attr_registry = defaultdict(make_tuple)

    def __getattr__(self, item):
        """
        Our method for Attrs is to register using the attribute type (i.e. Instrument) as keys in a dictionary.
        ``_attr_registry`` is a dictionary with the keys being subclasses of Attr
        and the value being the namedtuple of lists.
        As a result we index `_attr_registry` with `[self]` which will be the `type`
        of the `Attr` class to access the dictionary.
        This will return the namedtuple that has three attributes: `name`, `name_long` and `desc`.
        Each of which are a list.
        `name` will be the attribute name, `name_long` is the original name passed in and `desc` the description of the object.
        """
        # Get the revelant entries.
        registry = self._attr_registry[self]
        # All the attribute names under that type(Attr)
        names = registry.name
        if item in names:
            # We return Attr(name_long) to create the Attr requested.
            return self(registry.name_long[names.index(item)])
        else:
            raise AttributeError(f'This attribute, {item} is not defined, please register it.')

    def __dir__(self):
        """
        To tab complete in Python we need to add to the `__dir__()` return.
        So we add all the registered values for this subclass of Attr to the output.
        """
        return super().__dir__() + self._attr_registry[self].name

    def __str__(self):
        """
        This enables the "pretty" printing of Attrs.
        """
        name = ['Attribute Name'] + self._attr_registry[self].name
        name_long = ['Full Name'] + self._attr_registry[self].name_long
        desc = ['Description'] + self._attr_registry[self].desc
        pad_name = max(len(elm) for elm in name)
        pad_long = max(len(elm) for elm in name_long)
        pad_desc = max(len(elm) for elm in desc)
        fmt = "%-{}s | %-{}s | %-{}s".format(pad_name, pad_long, pad_desc)
        lines = [fmt % elm for elm in zip(name, name_long, desc)]
        lines.insert(1, '-' * (pad_name + 1) + '+' + '-' * (pad_long + 2) + '+' + '-' * (pad_desc + 1))
        return "\n".join(lines)


class Attr(metaclass=AttrMeta):
    """This is the base for all attributes."""
    def __and__(self, other):
        if isinstance(other, AttrOr):
            return AttrOr([elem & self for elem in other.attrs])
        if self.collides(other):
            return NotImplemented
        if isinstance(other, AttrAnd):
            return AttrAnd([self] + list(other.attrs))
        return AttrAnd([self, other])

    def __hash__(self):
        return hash(frozenset(vars(self).items()))

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
        Clients will use this method to register their values for subclasses of `~sunpy.net.attr.Attr`.

        The input has be a dictionary, with each key being a subclass of Attr.
        The value for each key should be a list of tuples with each tuple of the form (`Name`, `Description`).
        If you do not want to add a description, you can put `None` or an empty string.
        We sanitize the name you provide by removing all special characters and making it all lower case.
        If it still invalid we will append to the start of the name to make it a valid attribute name.
        This is an ugly workaround, so please try to pick a valid attribute name.

        **It is recommended that clients register values for all attrs used by their `_can_handle_query` method.**

        Parameters
        ----------

        adcit : `dict`
            A dictionary that has keys of `~sunpy.net.attr.Attr`.
            Each key should have a list of tuples as it value.
            Each tuple should be a pair of strings.
            First string is the attribute name you want to register
            Second string is the description of the attribute.

        Example
        -------
        >>> from sunpy.net import attr, attrs # doctest: +SKIP
        >>> attr.Attr.update_values({attrs.Instrument: [('AIA', 'AIA is in Space.'),
        ...                         ('HMI', 'HMI is next to AIA.')]}) # doctest: +SKIP
        >>> attr.Attr._attr_registry[attrs.Instrument] # doctest: +SKIP
        attr(name=['aia', 'hmi'], name_long=['AIA', 'HMI'],
            desc=['AIA is in Space.', 'HMI is next to AIA.'])
        >>> attr.Attr._attr_registry[attrs.Instrument].name # doctest: +SKIP
        ['aia', 'hmi']
        >>> attr.Attr._attr_registry[attrs.Instrument].name_long # doctest: +SKIP
        ['AIA', 'HMI']
        >>> attr.Attr._attr_registry[attrs.Instrument].desc # doctest: +SKIP
        ['AIA is in Space.', 'HMI is next to AIA.']
        >>> attrs.Instrument.aia, attrs.Instrument.hmi # doctest: +SKIP
        (<Instrument('AIA')>, <Instrument('HMI')>)
        """

        for key, value in adict.items():
            if isiterable(value) and not isinstance(value, str):
                for pair in value:
                    if len(pair) != 2:
                        raise ValueError(f'Invalid length (!=2) for values: {value}.')
                    else:
                        # Sanitize name, we remove all special characters and make it all lower case
                        name = ''.join(char for char in pair[0] if char.isalnum()).lower()
                        if keyword.iskeyword(name):
                            # Attribute name has been appended with `_`
                            # to make it a valid identifier since its a python keyword.
                            name = name + '_'
                        if not name.isidentifier():
                            # This should account for names with one number first.
                            # We match for single digits at the start only.
                            if _REGEX.match(name):
                                # This turns that digit into its name
                                number = strtonum(name[0])
                                name = number + ("_" + name[1:] if len(name) > 1 else "")
                        cls._attr_registry[key][0].append(name)
                        cls._attr_registry[key][1].append(pair[0])
                        cls._attr_registry[key][2].append(pair[1])
            else:
                raise ValueError((f"Invalid input value: {value} for key: {repr(key)}. "
                                  "The value is not iterable or just a string."))


class DummyAttr(Attr):
    """
    Empty attribute. Useful for building up queries. Returns other
    attribute when ORed or ANDed. It can be considered an empty query
    that you can use as an initial value if you want to build up your
    query in a loop.

    So, if we wanted an attr matching all the time intervals between the times
    stored as (from, to) tuples in a list, we could do.

    .. code-block:: python

       attr = DummyAttr()
       for from, to in times:
           attr |= Time(from, to)
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
        return hash(frozenset(self.attrs.items()))

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
