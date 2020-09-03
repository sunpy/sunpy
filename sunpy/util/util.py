"""
This module provides general utility functions.
"""
import os
import hashlib
from shutil import get_terminal_size
from itertools import chain, count
from collections import UserList

__all__ = ['unique', 'replacement_filename', 'expand_list',
           'expand_list_generator', 'dict_keys_same', 'hash_file', "get_width"]


def unique(itr, key=None):
    """
    Return only unique elements of a sequence.

    Parameters
    ----------
    itr : `iterable`
        Any iterable sequence.
    key : `function`, optional
        A function to apply to each element in the iterable. Defaults to `None`.

    Returns
    -------
    `set`:
        A `set` of each unique element.
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


def replacement_filename(path: str) -> str:
    """
    Return a replacement path if input path is currently in use.

    Enumerates until an unused filename is found, e.g., "foo.fits" becomes
    "foo.0.fits", if that is used, "foo.1.fits" and so on.

    Parameters
    ----------
    path : `str`
        A string path.

    Returns
    -------
    `str`:
        A string path.
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
    Expand a list of lists or tuples.

    Parameters
    ----------
    inp : `list`, `tuple`, `collections.UserList`
        The iterable to expand.

    Returns
    -------
    `list`
        A flat list consisting of the entries of the input.

    References
    ----------
    * https://stackoverflow.com/questions/2185822/expanding-elements-in-a-list/2185971#2185971
    """
    return [item for item in expand_list_generator(inp)]


def expand_list_generator(inp):
    for item in inp:
        if isinstance(item, (tuple, list, UserList)):
            yield from expand_list_generator(item)
        else:
            yield item


def partial_key_match(key, dictionary):
    """
    Return value/values from a dictionary based on a partial key.

    Each element of the partial key is matched against the keys of the dictionary and
    if a partial match is found the value of the key is returned.

    Even a partial match works here i.e even if the key matches partially a value is returned.

    Parameters
    ----------
    key : `tuple`
          A tuple containing the partial key.
    dictionary: `dict`
          The target dictionary from which we want to retrieve the value based on the partial key.

    Yields
    ------
    value
        The value of the matched key.

    References
    ----------
    * https://stackoverflow.com/questions/18893624/partial-match-dictionary-keyof-tuples-in-python

    Examples
    ----------
    >>> d = {('abc','def','ghi') : 1, ('abc', 'def', 'xyz') : 2, ('pqr', 'lmn', 'tuv') : 3}
    >>> list(partial_key_match(('abc', 'def', None), d))
        [1, 2]
    """
    for k, v in dictionary.items():
        if all(k1 == k2 or k2 is None for k1, k2 in zip(k, key)):
            yield v


def dict_keys_same(list_of_dicts):
    """
    Makes sure that a list of dictionaries all have the same keys.

    If a key is missing, it will be added but with a value of None.

    Parameters
    ----------
    list_of_dicts : `list` of `dict`
          A list containing each dictonary to parse.

    Returns
    ------
    `list`
        The list with each dict updated.

    References
    ----------
    * https://stackoverflow.com/questions/10482439/make-sure-all-dicts-in-a-list-have-the-same-keys

    Examples
    ----------
    >>> l = [{'x': 42}, {'x': 23, 'y': 5}]
    >>> dict_keys_same(l)
        [{'x': 42, 'y': None}, {'x': 23, 'y': 5}]
    """
    all_keys = set(chain.from_iterable(list_of_dicts))
    for d in list_of_dicts:
        d.update({key: None for key in all_keys if key not in d})
    return list_of_dicts


def hash_file(path):
    """
    Returns the SHA-256 hash of a file.

    Parameters
    ----------
    path: `str`
        The path of the file to be hashed.

    Returns
    -------
    `str`
        SHA-256 hash of the file.

    References
    ----------
    * https://stackoverflow.com/a/22058673
    """
    BUF_SIZE = 65536
    sha256 = hashlib.sha256()

    with open(path, 'rb') as f:
        while True:
            data = f.read(BUF_SIZE)
            if not data:
                break
            sha256.update(data)

    return sha256.hexdigest()


def get_width():
    """
    Gets the width of the current terminal.
    Accounts for if the 'COLUMNS' environmental variable is set.

    Returns
    -------
    `int`
        Width of the terminal you are in.
        Works for IPython notebooks and normal terminals.
    """
    width = os.environ.get("COLUMNS", None)
    if width:
        width = int(width)
    else:
        width, _ = get_terminal_size()
    return width
