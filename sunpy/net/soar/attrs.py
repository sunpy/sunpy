import warnings

import sunpy.net.attrs as a
from sunpy.net.attr import AttrAnd, AttrOr, AttrWalker, DataAttr, SimpleAttr
from sunpy.util.exceptions import SunpyDeprecationWarning

__all__ = ['Product']


class Product(SimpleAttr):
    """
    The data product identifier to search for.
    """


class Identifier(Product):
    """
    The data product identifier to search for.
    """
    def __init__(self, *args, **kwargs):
        warnings.warn(
            "'a.soar.Identifier' is deprecated; use 'a.soar.Product' instead.",
            SunpyDeprecationWarning)
        super().__init__(*args, **kwargs)


walker = AttrWalker()


@walker.add_creator(AttrOr)
def create_or(wlk, tree):
    """
    Creator for OR. Loops through the next level down in the tree and appends
    the individual results to a list.
    """
    results = []
    for sub in tree.attrs:
        results.append(wlk.create(sub))
    return results


@walker.add_creator(AttrAnd, DataAttr)
def create_and(wlk, tree):
    """
    Creator for And and other simple attributes. No walking needs to be done,
    so simply call the applier function.
    """
    result = []
    wlk.apply(tree, result)
    return [result]


@walker.add_applier(AttrAnd)
def apply_and(wlk, and_attr, params):
    """
    Applier for And.

    Parameters
    ----------
    wlk : AttrWalker
    and_attr : AttrAnd
        The AND attribute being applied. The individual attributes being
        AND'ed together are accesible with ``and_attr.attrs``.
    params : list[str]
        List of search parameters.
    """
    for iattr in and_attr.attrs:
        wlk.apply(iattr, params)


"""
Below are appliers for individual attributes.

The all convert the attribute object into a query string, that will eventually
be passed as a query to the SOAR server. They all have the signature:

Parameters
----------
wlk : AttrWalker
    The attribute walker.
attr :
    The attribute being applied.
params : list[str]
    List of search parameters.
"""


@walker.add_applier(a.Time)
def _(wlk, attr, params):
    start = attr.start.strftime('%Y-%m-%d+%H:%M:%S')
    end = attr.end.strftime('%Y-%m-%d+%H:%M:%S')
    params.append(f"begin_time>='{start}'+AND+begin_time<='{end}'")


@walker.add_applier(a.Level)
def _(wlk, attr, params):
    level = int(attr.value)
    params.append(f"level='L{level}'")


@walker.add_applier(a.Instrument)
def _(wlk, attr, params):
    params.append(f"instrument='{attr.value}'")


@walker.add_applier(Product, Identifier)
def _(wlk, attr, params):
    params.append(f"descriptor='{attr.value}'")
