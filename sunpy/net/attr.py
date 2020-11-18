"""
Allow representation of queries as logic expressions. This module makes
sure that attributes that are combined using the two logic operations AND (&)
and OR (|) always are in disjunctive normal form, that is, there are only two
levels ­- the first being disjunction and the second being conjunction. In other
words, every combinations of attributes looks like this:
(a AND b AND c) OR (d AND e).

Walkers are used to traverse the tree that results from combining attributes.
They are implemented using `functools.singledispatch` modified to dispatch on the second argument to the function.

Please note that & is evaluated first, so A & B | C is equivalent to
(A & B) | C.
"""
import re
import string
import keyword
import textwrap
from textwrap import dedent
from collections import namedtuple, defaultdict

from astropy.table import Table
from astropy.utils.misc import isiterable

from sunpy.extern import inflect
from sunpy.util.functools import seconddispatch
from sunpy.util.util import get_width

_ATTR_TUPLE = namedtuple("attr", "name client name_long desc")
# Matches any number.
NUMBER_REGEX = re.compile(r"^(\d+$|\d(?:\.\d+)?)")

__all__ = ['Attr', 'DummyAttr', 'SimpleAttr', 'Range', 'AttrAnd', 'AttrOr',
           'ValueAttr', 'and_', 'or_', 'AttrWalker']


def make_tuple():
    return _ATTR_TUPLE([], [], [], [])


def _print_attrs(attr, html=False):
    """
    Given a Attr class will print out each registered attribute.

    Parameters
    ----------
    attr : `sunpy.net.attr.Attr`
        The attr class/type to print for.
    html : bool
        Will return a html table instead.

    Returns
    -------
    `str`
        String with the registered attributes.
    """
    attrs = attr._attr_registry[attr]
    # Only sort the attrs if any have been registered
    sorted_attrs = _ATTR_TUPLE(*zip(*sorted(zip(*attrs)))) if attrs.name else make_tuple()
    *other_row_data, descs = sorted_attrs
    descs = [(dsc[:77] + '...') if len(dsc) > 80 else dsc for dsc in descs]
    table = Table(names=["Attribute Name", "Client", "Full Name", "Description"],
                  dtype=["U80", "U80", "U80", "U80"],
                  data=[*other_row_data, descs])

    class_name = f"{(attr.__module__ + '.') or ''}{attr.__name__}"
    lines = [class_name]
    # If the attr lacks a __doc__ this will error and prevent this from returning anything.
    try:
        lines.append(dedent(attr.__doc__.partition("\n\n")[0]) + "\n")
    except AttributeError:
        pass

    format_line = "<p>{}</p>" if html else "{}"
    width = -1 if html else get_width()

    lines = [*[format_line.format(line) for line in lines],
             *table.pformat_all(show_dtype=False, max_width=width, align="<", html=html)]
    return '\n'.join(lines)


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
        try:
            # We return Attr(name_long) to create the Attr requested.
            return self(registry.name_long[names.index(item)])
        except ValueError:
            raise AttributeError(f'This attribute, {item} is not defined, please register it.')

    def __dir__(self):
        """
        To tab complete in Python we need to add to the `__dir__()` return.
        So we add all the registered values for this subclass of Attr to the output.
        """
        custom_attrs = set(self._attr_registry[self].name)
        # "all" can be registered as a documentation helper, but isn't a valid attr
        custom_attrs.discard("all")
        return super().__dir__() + list(custom_attrs)

    def __repr__(self):
        """
        Returns the normal repr plus the pretty attr __str__.
        """
        return f"{type.__repr__(self)}\n{str(self)}"

    def __str__(self):
        """
        This enables the "pretty" printing of Attrs.
        """
        return _print_attrs(self)

    def _repr_html_(self):
        """
        This enables the "pretty" printing of Attrs with html.
        """
        return _print_attrs(self, html=True)


class Attr(metaclass=AttrMeta):
    """This is the base for all attributes."""

    def __and__(self, other):
        if isinstance(other, AttrOr):
            return AttrOr([elem & self for elem in other.attrs])
        if self.collides(other):
            return NotImplemented
        if isinstance(other, AttrAnd):
            return AttrAnd([self, *other.attrs])
        return AttrAnd([self, other])

    def __hash__(self):
        return hash(frozenset(vars(self).items()))

    def __or__(self, other):
        # Optimization.
        if self == other:
            return self
        if isinstance(other, AttrOr):
            return AttrOr([self, *other.attrs])
        return AttrOr([self, other])

    def collides(self, other):
        raise NotImplementedError

    def __eq__(self, other):
        if not isinstance(other, Attr):
            return False
        return vars(self) == vars(other)

    @classmethod
    def update_values(cls, adict):
        """
        Clients will use this method to register their values for subclasses of `~sunpy.net.attr.Attr`.

        The input has to be a dictionary, with each key being an instance of a client.
        The value for each client has to be a dictionary with each key being a subclass of Attr.
        The value for each Attr key should be a list of tuples with each tuple of the form (`Name`, `Description`).
        If you do not want to add a description, you can put `None` or an empty string.
        We sanitize the name you provide by removing all special characters and making it all lower case.
        If it still invalid we will append to the start of the name to make it a valid attribute name.
        This is an ugly workaround, so please try to pick a valid attribute name.

        **It is recommended that clients register values for all attrs used by their `_can_handle_query` method.**

        Parameters
        ----------
        adict : `dict`
            The keys for this dictionary have to be an instance of a downloader client.
            The values for each client are a dictionary that has keys of `~sunpy.net.attr.Attr`.
            Each key should have a list of tuples as it value.
            Each tuple should be a pair of strings.
            First string is the attribute name you want to register
            Second string is the description of the attribute.

        Example
        -------
        # The first import is to make this example work, it should not be used otherwise
        >>> from sunpy.net.dataretriever import GenericClient
        >>> from sunpy.net import attr, attrs
        >>> attr.Attr.update_values({GenericClient : {
        ...                          attrs.Instrument: [('AIA', 'AIA is in Space.'),
        ...                                             ('HMI', 'HMI is next to AIA.')]}})   # doctest: +SKIP
        >>> attr.Attr._attr_registry[attrs.Instrument]  # doctest: +SKIP
        attr(name=['aia', 'hmi'], client=['Generic', 'Generic'], name_long=['AIA', 'HMI'],
        desc=['AIA is in Space.', 'HMI is next to AIA.'])
        >>> attr.Attr._attr_registry[attrs.Instrument].name  # doctest: +SKIP
        ['aia', 'hmi']
        >>> attr.Attr._attr_registry[attrs.Instrument].name_long  # doctest: +SKIP
        ['AIA', 'HMI']
        >>> attr.Attr._attr_registry[attrs.Instrument].desc  # doctest: +SKIP
        ['AIA is in Space.', 'HMI is next to AIA.']
        >>> attrs.Instrument.aia, attrs.Instrument.hmi  # doctest: +SKIP
        (<sunpy.net.attrs.Instrument(AIA: AIA is in Space.) object at 0x...>,
        <sunpy.net.attrs.Instrument(HMI: HMI is next to AIA.) object at 0x...>)
        """
        for client, attr_dict in adict.items():
            for attr, attr_values in attr_dict.items():
                if not isiterable(attr_values) or isinstance(attr_values, str):
                    raise ValueError(f"Invalid input value: {attr_values} for key: {repr(attr)}. "
                                     "The value is not iterable or just a string.")

                attr_tuple = cls._attr_registry[attr]

                for pair in attr_values:
                    if len(pair) > 2:
                        raise ValueError(f'Invalid length (!=2) for values: {attr_values}.')
                    elif len(pair) == 1:
                        if pair[0] != "*":
                            raise ValueError(
                                f'Invalid value given for * registration: {attr_values}.')
                        # Special case handling for * aka all values allowed.
                        pair = ["all", "All values of this type are supported."]

                    # Sanitize part one: Check if the name has a number in it
                    number_match = NUMBER_REGEX.match(pair[0])
                    p = inflect.engine()
                    try:
                        number_str = number_match.group(1)
                        name = p.number_to_words(number_str)
                        if number_str != number_match.string:
                            name = name + "_" + number_match.string[number_match.end(1):]
                    except AttributeError:
                        name = pair[0]

                    # Sanitize part two: remove punctuation and replace it with _
                    name = re.sub('[%s]' % re.escape(string.punctuation), '_', name)
                    # Sanitize name, we remove all special characters
                    name = ''.join(char for char in name
                                   if char.isidentifier() or char.isnumeric())
                    # Make name lower case
                    name = name.lower()

                    if keyword.iskeyword(name):
                        # Attribute name has been appended with `_`
                        # to make it a valid identifier since its a python keyword.
                        name = name + '_'

                    attr_tuple[0].append(name)
                    attr_tuple[1].append(client.__name__.replace("Client", ""))
                    attr_tuple[2].append(pair[0])
                    attr_tuple[3].append(pair[1])


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
        obj_placeholder = " object "
        attr_reg = AttrMeta._attr_registry[self.__class__]
        new_repr = object.__repr__(self).split(obj_placeholder)
        # If somehow the idx isn't in the attr reg,
        # we still want it to print it repr without error.
        try:
            idx = attr_reg.name_long.index(self.value)
            obj_value_repr = f"({self.value}: {attr_reg.desc[idx]})"
        except ValueError:
            obj_value_repr = f": {self.value}"
        new_repr = [new_repr[0], obj_value_repr, obj_placeholder, new_repr[1]]
        return textwrap.fill("".join(new_repr), 100)

    @property
    def type_name(self):
        return self.__class__.__name__.lower()


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
            return AttrAnd([*self.attrs, *other.attrs])
        if isinstance(other, AttrOr):
            return AttrOr([elem & self for elem in other.attrs])
        return AttrAnd([*self.attrs, other])

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
            return AttrOr([*self.attrs, *other.attrs])
        return AttrOr([*self.attrs, other])

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

    def __repr__(self):
        creators = list(self.createmm.registry.keys())
        appliers = list(self.applymm.registry.keys())
        return f"""{super().__repr__()}
Registered creators: {creators}
Registered appliers: {appliers}"""

    @staticmethod
    def _unknown_type(*args, **kwargs):
        raise TypeError(
            f"{args[1]} or any of its parents have not been registered with the AttrWalker")

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
