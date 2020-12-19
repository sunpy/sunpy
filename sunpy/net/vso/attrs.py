# Author: Florian Mayer <florian.mayer@bitsrc.org>
#
# This module was developed with funding provided by
# the ESA Summer of Code (2011).
#
"""
Attributes that can be used to construct VSO queries.

Attributes are the fundamental building blocks of queries that, together with
the two operations of AND and OR (and in some rare cases XOR) can be used to
construct complex queries. Most attributes can only be used once in an
AND-expression, if you still attempt to do so it is called a collision. For a
quick example think about how the system should handle Instrument('aia') &
Instrument('eit').
"""


import sys
import warnings
from functools import singledispatch as _singledispatch

import astropy.units as u
from astropy.time import Time as AstropyTime

from sunpy.util.exceptions import SunpyDeprecationWarning
from .. import _attrs
from .. import attr as _attr

__all__ = ['Extent', 'Field', 'Pixels', 'Filter', 'Quicklook', 'PScale']


# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__


_TIMEFORMAT = '%Y%m%d%H%M%S'


class Field(_attr.ValueAttr):
    """
    A subclass of the value attribute.  Used in defining a decorator for the
    dummy attribute.
    """

    def __init__(self, fielditem):
        _attr.ValueAttr.__init__(self, {
            ('field', 'fielditem'): fielditem
        })


class Extent(_attr.DataAttr):
    """
    Specify the spatial field-of-view of the query. Due to a bug in the VSO,
    the Extent attribute is not used.

    """

    def __init__(self, x, y, width, length, atype):
        super().__init__()

        self.x = x
        self.y = y
        self.width = width
        self.length = length
        self.type = atype

    def collides(self, other):
        return isinstance(other, self.__class__)


class Pixels(_attr.SimpleAttr):
    """
    Pixels are (currently) limited to a single dimension (and only implemented
    for SDO data)  We hope to change this in the future to support TRACE,
    Hinode and other investigations where this changed between observations.

    References
    ----------
    Documentation in SSWIDL routine vso_search.pro.
    """


class PScale(_attr.SimpleAttr):
    """
    Pixel Scale (PSCALE) is in arc seconds.

    Parameters
    ----------
    value : float or str

        The value can be entered in of three ways:

        #. May be entered as a string or any numeric type for equality matching
        #. May be a string of the format '(min) - (max)' for range matching
        #. May be a string of the form '(operator) (number)' where operator is\
        one of: lt gt le ge < > <= >=

        Currently only implemented for SDO, which is 0.6 arcsec per pixel at full
        resolution for AIA.

    References
    ----------
    Documentation in SSWIDL routine vso_search.pro.
    """


class Quicklook(_attr.SimpleAttr):
    """
    Retrieve 'quicklook' data if available.

    Parameters
    ----------
    value : boolean
        Set to True to retrieve quicklook data if available.

        Quicklook items are assumed to be generated with a focus on speed rather
        than scientific accuracy.  They are useful for instrument planning and
        space weather but should not be used for science publication.
        This concept is sometimes called 'browse' or 'near real time' (nrt)
        Quicklook products are *not* searched by default.

    References
    ----------
    Documentation in SSWIDL routine vso_search.pro.
    """

    def __init__(self, value):
        super().__init__(value)
        if self.value:
            self.value = 1
        else:
            self.value = 0


class Filter(_attr.SimpleAttr):
    """
    This attribute is a placeholder for the future.

    Parameters
    ----------
    value : str

    """


# The walker specifies how the Attr-tree is converted to a query the
# server can handle.
_walker = _attr.AttrWalker()

# The _create functions make a new VSO query from the attribute tree,
# the _apply functions take an existing query-block and update it according
# to the attribute tree passed in as root. Different attributes require
# different functions for conversion into query blocks.


@_walker.add_creator(_attr.ValueAttr, _attr.AttrAnd)
def _create(wlk, root, api):
    """ Implementation detail. """
    api.set_ns_prefix('VSO', 'http://virtualsolar.org/VSO/VSOi')
    value = api.get_type('VSO:QueryRequestBlock')()
    wlk.apply(root, api, value)
    return [value]


@_walker.add_applier(_attr.ValueAttr)
def _apply(wlk, root, api, block):
    """ Implementation detail. """
    for k, v in root.attrs.items():
        name = k[0]
        subkey = k[1:]

        if subkey:
            if len(subkey) != 1:
                raise ValueError("Can't parse double nested ValueAttr")
            subkey = subkey[0]

            if block[name]:
                block[name].update({subkey: v})
            else:
                block[name] = {subkey: v}
        else:
            block[name] = v


@_walker.add_applier(_attr.AttrAnd)
def _apply(wlk, root, api, queryblock):
    """ Implementation detail. """
    for attr in root.attrs:
        wlk.apply(attr, api, queryblock)


@_walker.add_creator(_attr.AttrOr)
def _create(wlk, root, api):
    """ Implementation detail. """
    blocks = []
    for attr in root.attrs:
        blocks.extend(wlk.create(attr, api))
    return blocks


# Converters take a type unknown to the walker and convert it into one
# known to it. All of those convert types into ValueAttrs, which are
# handled above by just assigning according to the keys and values of the
# attrs member.
_walker.add_converter(Extent)(
    lambda x: _attr.ValueAttr(
        {('extent', k): v for k, v in vars(x).items()}
    )
)

_walker.add_converter(_attrs.Time)(
    lambda x: _attr.ValueAttr({
        ('time', 'start'): x.start.strftime(_TIMEFORMAT),
        ('time', 'end'): x.end.strftime(_TIMEFORMAT),
        ('time', 'near'): (
            x.near.strftime(_TIMEFORMAT) if x.near is not None else None),
    })
)

_walker.add_converter(_attr.SimpleAttr)(
    lambda x: _attr.ValueAttr({(x.__class__.__name__.lower(), ): x.value})
)

_walker.add_converter(_attrs.Wavelength)(
    lambda x: _attr.ValueAttr({
        ('wave', 'wavemin'): x.min.value,
        ('wave', 'wavemax'): x.max.value,
        ('wave', 'waveunit'): x.unit.name,
    })
)

# The idea of using a multi-method here - that means a method which dispatches
# by type but is not attached to said class - is that the attribute classes are
# designed to be used not only in the context of VSO but also elsewhere (which
# AttrAnd and AttrOr obviously are - in the HEK module). If we defined the
# filter method as a member of the attribute classes, we could only filter
# one type of data (that is, VSO data).


@_singledispatch
def _filter_results(*args, **kwargs):
    raise TypeError(f"a[0] is not a type that can filter QueryResponse objects.")


# If we filter with ANDed together attributes, the only items are the ones
# that match all of them - this is implementing  by ANDing the pool of items
# with the matched items - only the ones that match everything are there
# after this.
@_filter_results.register(_attr.AttrAnd)
def _(attr, results):
    res = set(results)
    for elem in attr.attrs:
        res &= _filter_results(elem, res)
    return res


# If we filter with ORed attributes, the only attributes that should be
# removed are the ones that match none of them. That's why we build up the
# resulting set by ORing all the matching items.
@_filter_results.register(_attr.AttrOr)
def _(attr, results):
    res = set()
    for elem in attr.attrs:
        res |= _filter_results(elem, results)
    return res


# Filter out items by comparing attributes.
@_filter_results.register(_attr.SimpleAttr)
def _(attr, results):
    attrname = attr.__class__.__name__.lower()
    return {
        item for item in results
        # Some servers seem to omit some fields. No way to filter there.
        if not hasattr(item, attrname) or
        getattr(item, attrname).lower() == attr.value.lower()
    }


# The dummy attribute does not filter at all.
@_filter_results.register(_attr.DummyAttr, Field)
def _(attr, results):
    return set(results)


@_filter_results.register(_attrs.Wavelength)
def _(attr, results):
    return {
        it for it in results
        if
        it.wave.wavemax is not None
        and
        attr.min <= u.Quantity(it.wave.wavemax, unit=it.wave.waveunit).to(
            u.angstrom, equivalencies=u.spectral())
        and
        it.wave.wavemin is not None
        and
        attr.max >= u.Quantity(it.wave.wavemin, unit=it.wave.waveunit).to(
            u.angstrom, equivalencies=u.spectral())
    }


@_filter_results.register(_attrs.Time)
def _(attr, results):
    return {
        it for it in results
        if
        it.time.end is not None
        and
        attr.min <= AstropyTime.strptime(it.time.end, _TIMEFORMAT)
        and
        it.time.start is not None
        and
        attr.max >= AstropyTime.strptime(it.time.start, _TIMEFORMAT)
    }


@_filter_results.register(Extent)
def _(attr, results):
    return {
        it for it in results
        if
        it.extent.type is not None
        and
        it.extent.type.lower() == attr.type.lower()
    }


# Deprecate old classes
class _DeprecatedAttr:
    def __init__(self, *args, **kwargs):
        name = type(self).__name__
        warnings.warn(f"sunpy.net.vso.attrs.{name} is deprecated, please use sunpy.net.attrs.{name}",
                      SunpyDeprecationWarning)
        super().__init__(*args, **kwargs)


_deprecated_names = ['Time', 'Instrument', 'Wavelength', 'Source', 'Provider',
                     'Level', 'Sample', 'Detector', 'Resolution', 'Physobs']

for _name in _deprecated_names:
    # Dynamically construct a class which inherits the class with the
    # deprecation warning in the __init__ first and the class it's deprecating
    # second.
    cls = type(_name, (_DeprecatedAttr, getattr(_attrs, _name)), {})
    # Add the new class to the modules namespace
    setattr(sys.modules[__name__], _name, cls)

__all__ += _deprecated_names
