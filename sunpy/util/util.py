"""
General utility functions.
"""

from __future__ import absolute_import, division, print_function

import os
from itertools import count

import numpy as np

from sunpy.extern import six
from sunpy.extern.six.moves import map, zip

__all__ = ['to_signed', 'unique', 'print_table', 'replacement_filename',
           'merge', 'common_base', 'minimal_pairs', 'expand_list',
           'expand_list_generator']


def to_signed(dtype):
    """
    Return dtype that can hold data of passed dtype but is signed.
    Raise ValueError if no such dtype exists.

    Parameters
    ----------
    dtype : `numpy.dtype`
        dtype whose values the new dtype needs to be able to represent.

    Returns
    -------
    `numpy.dtype`

    """
    if dtype.kind == "u":
        if dtype.itemsize == 8:
            raise ValueError("Cannot losslessly convert uint64 to int.")
        dtype = "int{0:d}".format(min(dtype.itemsize * 2 * 8, 64))
    return np.dtype(dtype)


def unique(itr, key=None):
    """
    not documented yet

    Parameters
    ----------
    itr : iterable
        Object to be iterated over

    key : object
        not documented yet
    """
    items = set()
    if key is None:
        for elem in itr:
            if elem not in items:
                yield elem
                items.add(elem)
    else:
        for elem in itr:
            x = key(elem)
            if x not in items:
                yield elem
                items.add(x)


def print_table(lst, colsep=' ', linesep='\n'):
    width = [max(map(len, col)) for col in zip(*lst)]
    return linesep.join(
        colsep.join(col.ljust(n) for n, col in zip(width, row)) for row in lst)


def minimal_pairs(one, other):
    """ Find pairs of values in one and other with minimal distance.
    Assumes one and other are sorted in the same sort sequence.

    Parameters
    ----------
    one, other : sequence
        Sequence of scalars to find pairs from.

    Returns
    -------
    `tuple`
         Pairs of values in `one` and `other` with minimal distance
    """
    lbestdiff = bestdiff = bestj = besti = None
    for i, freq in enumerate(one):
        lbestj = bestj

        bestdiff, bestj = None, None
        for j, o_freq in enumerate(other[lbestj:]):
            j = lbestj + j if lbestj else j
            diff = abs(freq - o_freq)
            if bestj is not None and diff > bestdiff:
                break

            if bestj is None or bestdiff > diff:
                bestj = j
                bestdiff = diff

        if lbestj is not None and lbestj != bestj:
            yield (besti, lbestj, lbestdiff)
            besti = i
            lbestdiff = bestdiff
        elif lbestdiff is None or bestdiff < lbestdiff:
            besti = i
            lbestdiff = bestdiff

    yield (besti, bestj, lbestdiff)


DONT = object()


def find_next(one, other, pad=DONT):
    """
    Given two sorted sequences one and other, for every element
    in one, return the one larger than it but nearest to it in other.
    If no such exists and pad is not DONT, return value of pad as "partner".
    """
    n = 0
    for elem1 in one:
        for elem2 in other[n:]:
            n += 1
            if elem2 > elem1:
                yield elem1, elem2
                break
        else:
            if pad is not DONT:
                yield elem1, pad


def common_base(objs):
    """
    Find class that every item of objs is an instance of.
    """
    for cls in objs[0].__class__.__mro__:
        if all(isinstance(obj, cls) for obj in objs):
            break
    return cls


def merge(items, key=(lambda x: x)):
    """
    Given sorted lists of iterables, return new iterable that returns
    elements of all iterables sorted with respect to key.
    """
    state = {}
    for item in map(iter, items):
        try:
            first = next(item)
        except StopIteration:
            continue
        else:
            state[item] = (first, key(first))

    while state:
        for item, (value, tk) in six.iteritems(state):
            # Value is biggest.
            if all(tk >= k for it, (v, k) in six.iteritems(state)
                   if it is not item):
                yield value
                break
        try:
            n = next(item)
            state[item] = (n, key(n))
        except StopIteration:
            del state[item]


def replacement_filename(path):
    """
    Return replacement path for already used path. Enumerates
    until an unused filename is found. E.g., "/home/florian/foo.fits"
    becomes "/home/florian/foo.0.fits", if that is used
    "/home/florian/foo.1.fits", etc.
    """
    if not os.path.exists(path):
        return path
    else:
        dir_, filename = os.path.split(path)
        base, ext = os.path.splitext(filename)
        for c in count():
            name = base + '.' + str(c) + ext
            newpath = os.path.join(dir_, name)
            if not os.path.exists(newpath):
                return newpath


def expand_list(inp):
    """
    Expand a list of lists.

    Parameters
    ----------
    inp : `list`

    Returns
    -------
    `list`
        A flat list consisting of the entries of the input.

    References
    ----------
    Taken from :https://stackoverflow.com/questions/2185822/expanding-elements-in-a-list/2185971#2185971

    """
    return [item for item in expand_list_generator(inp)]


def expand_list_generator(inp):
    for item in inp:
        if type(item) in [list, tuple]:
            for nested_item in expand_list_generator(item):
                yield nested_item
        else:
            yield item
