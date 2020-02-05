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

from sunpy.util.functools import seconddispatch

_ATTR_TUPLE = namedtuple("attr", "name name_long desc")
_REGEX = re.compile(r"^[\d]([^\d].*)?$")

__all__ = ['Attr', 'DummyAttr', 'SimpleAttr', 'Range', 'AttrAnd', 'AttrOr',
           'ValueAttr', 'and_', 'or_', 'AttrWalker']


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
        fmt = f"%-{pad_name}s | %-{pad_long}s | %-{pad_desc}s"
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
                raise ValueError(f"Invalid input value: {value} for key: {repr(key)}. "
                                  "The value is not iterable or just a string.")


class DataAttr(Attr):
    """
    A base class for attributes classes which contain data.

    This is to differentiate them from classes like `sunpy.net.attr.AttrAnd` or
    the base `sunpy.net.attr.Attr` class which do not. The motivation for this
    distinction is to make it easier for walkers to match all classes which are
    not user specified Attrs.
    """
    def __new__(cls, *args, **kwargs):
        if cls is DataAttr:
            raise TypeError("You should not directly instantiate DataAttr, only it's subclasses.")

        return super().__new__(cls)


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


class SimpleAttr(DataAttr):
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
        super().__init__()
        self.value = value

    def collides(self, other):
        return isinstance(other, self.__class__)

    def __repr__(self):
        return "<{cname!s}({val!r})>".format(
            cname=self.__class__.__name__, val=self.value)


class Range(DataAttr):
    """
    An attribute that represents a range of a value.

    This type of attribute would be applicable for types like Wavelength or Time.
    The range is inclusive of both the min and ma

    Parameters
    ----------
    min_ : `object`
        The lower bound of the range.
    max_ : `object`
        The upper bound of the range.
    """
    def __init__(self, min_, max_):
        self.min = min_
        self.max = max_

        super().__init__()

    def __xor__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented

        new = DummyAttr()
        if self.min < other.min:
            new |= type(self)(self.min, min(other.min, self.max))
        if other.max < self.max:
            new |= type(self)(other.max, self.max)
        return new

    def __contains__(self, other):
        if isinstance(other, Range):
            return self.min <= other.min and self.max >= other.max
        else:
            return self.min <= other <= self.max


class AttrAnd(Attr):
    """ Attribute representing attributes ANDed together. """
    def __init__(self, attrs):
        super().__init__()
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
        return f"<AttrAnd({self.attrs!r})>"

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
        super().__init__()
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
        return f"<AttrOr({self.attrs!r})>"

    def __eq__(self, other):
        if not isinstance(other, AttrOr):
            return False
        return set(self.attrs) == set(other.attrs)

    def __hash__(self):
        return hash(frozenset(self.attrs))

    def collides(self, other):
        return all(elem.collides(other) for elem in self.attrs)


# This appears to only be used as a base type for the Walker, i.e. a common
# denominator for the walker to convert to whatever the output of the walker is
# going to be.
class ValueAttr(DataAttr):
    def __init__(self, attrs):
        super().__init__()
        self.attrs = attrs

    def __repr__(self):
        return f"<ValueAttr({self.attrs!r})>"

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


class AttrWalker:
    """
    Traverse the Attr tree and convert it to a different representation.

    The ``AttrWalker`` can walk a complex tree of attrs and represent that tree
    in a way that is useful to the client using the attrs. For the VSO client
    it generates a ``VSO:QueryResponseBlock`` object, for the database module
    it performs database queries and returns results from the database.

    The walker has three core operations that can be applied to the tree, all
    of these are functions which are applied to one or more
    `~sunpy.net.attr.Attr` types, using conditional dispatch based on type.

    * creators: Creators when given an `~sunpy.net.attr.Attr` return a new object.
    * appliers: Appliers process an `~sunpy.net.attr.Attr` type and modify any
      arguments passed.
    * converters: Converters convert types unknown to any other creator or
      appliers to types known by them. They take in an `~sunpy.net.attr.Attr`
      type and return a different one.

    """
    @staticmethod
    def _unknown_type(*args, **kwargs):
        raise TypeError(f"{args[1]} or any of its parents have not been registered with the AttrWalker")

    def __init__(self):
        self.applymm = seconddispatch(self._unknown_type)
        self.createmm = seconddispatch(self._unknown_type)

    def create(self, *args, **kwargs):
        """
        Call the create function(s) matching the arguments to this method.
        """
        return self.createmm(self, *args, **kwargs)

    def apply(self, *args, **kwargs):
        """
        Call the apply function(s) matching the arguments to this method.
        """
        return self.applymm(self, *args, **kwargs)

    def add_creator(self, *types):
        """
        Register all specified types with this function for the ``.create`` method.
        """
        def _dec(fun):
            for type_ in types:
                self.createmm.register(type_, fun)
            return fun
        return _dec

    def add_applier(self, *types):
        """
        Register all specified types with this function for the ``.apply`` method.
        """
        def _dec(fun):
            for type_ in types:
                self.applymm.register(type_, fun)
            return fun
        return _dec

    def add_converter(self, *types):
        """
        Register a function to convert the specified type into a known type for create and apply.

        After a converter is run, create or apply will be called again with the new types.
        """
        def _dec(fun):
            for type_ in types:
                self.applymm.register(type_, self._cv_apply(fun))
                self.createmm.register(type_, self._cv_create(fun))
            return fun
        return _dec

    def _cv_apply(self, fun):
        """
        Call the converter and then re-call apply.
        """
        def _fun(*args, **kwargs):
            args = list(args)
            args[1] = fun(args[1])
            return self.applymm(*args, **kwargs)
        return _fun

    def _cv_create(self, fun):
        """
        Call the converter and then re-call convert.
        """
        def _fun(*args, **kwargs):
            args = list(args)
            args[1] = fun(args[1])
            return self.createmm(*args, **kwargs)
        return _fun


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
