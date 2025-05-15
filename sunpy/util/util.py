"""
This module provides general utility functions.
"""
import os
import hashlib
import inspect
from io import BytesIO
from base64 import b64encode
from shutil import get_terminal_size
from itertools import chain, count
from collections import UserList
from collections.abc import Iterator

import numpy as np

from astropy.wcs.utils import pixel_to_pixel

__all__ = ['unique', 'replacement_filename', 'expand_list',
           'expand_list_generator', 'dict_keys_same', 'hash_file', 'get_width',
           'get_keywords', 'get_set_methods', 'fix_duplicate_notes',
           'grid_perimeter', 'extent_in_other_wcs']

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


def replacement_filename(path: str):
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
    return list(expand_list_generator(inp))


def expand_list_generator(inp):
    for item in inp:
        # parfive.Results are UserList
        if isinstance(item, (list | tuple | UserList | Iterator)):
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
    dictionary : `dict`
        The target dictionary from which we want to retrieve the value based on the partial key.

    Yields
    ------
    value
        The value of the matched key.

    References
    ----------
    * https://stackoverflow.com/questions/18893624/partial-match-dictionary-keyof-tuples-in-python

    Examples
    --------
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
          A list containing each dictionary to parse.

    Returns
    ------
    `list`
        The list with each dict updated.

    References
    ----------
    * https://stackoverflow.com/questions/10482439/make-sure-all-dicts-in-a-list-have-the-same-keys

    Examples
    --------
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
    path : `str`
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


def get_keywords(func):
    """
    Returns a set of keyword names from ``func``'s signature.
    Recursive if ``func`` is a list of functions and methods.

    Parameters
    ----------
    func : function or method or `list`
        Function or method (or list of those) to extract a
        set of accepted keyword arguments for.

    Returns
    -------
    keywords : `set`
        A set of accepted keyword arguments.
    """
    if isinstance(func, list):
        keywords = set()
        for f in func:
            keywords.update(get_keywords(f))
        return keywords
    sig = inspect.signature(func)
    keywords = {param.name for param in sig.parameters.values()
                if param.default is not inspect.Parameter.empty}
    return keywords


def get_set_methods(obj):
    """
    Returns a set of keyword names that can be handled by
    an object's ``set_...`` methods.

    Parameters
    ----------
    obj : `object`
        Matplotlib object such as `~matplotlib.image.AxesImage`
        to extract handled keyword arguments for.

    Returns
    -------
    keywords : `set`
        A set of accepted keyword arguments.

    Notes
    -----
    See :meth:`matplotlib.artist.Artist.update` for an example
    of Matplotlib relying on this capability.
    """
    return {
        m[4:] for m in dir(obj)
        if m.startswith('set_') and callable(getattr(obj, m, None))
    }


def _figure_to_base64(fig):
    # Converts a matplotlib Figure to a base64 UTF-8 string
    buf = BytesIO()
    fig.savefig(buf, format='png', facecolor='none')  # works better than transparent=True
    return b64encode(buf.getvalue()).decode('utf-8')


def fix_duplicate_notes(notes_to_add, docstring):
    """
    Merges a note section into a docstring that may have a notes section.

    Parameters
    ----------
    notes_to_add : str
        The notes that need to be added (starting with ``"Notes"``).
    docstring
        The original docstring.

    Returns
    -------
    str
        Updated docstring with a single notes section.
    """
    docstring = inspect.cleandoc(docstring)  # not necessary on Python 3.13+
    existing_notes_pos = docstring.find('Notes\n-----')
    new_notes_pos = notes_to_add.find('Notes\n-----')
    if new_notes_pos == -1:
        raise RuntimeError("The notes to add does not appear to be a 'Notes' section.")
    notes_to_add_data = notes_to_add[new_notes_pos + len('Notes\n-----'):].strip()

    # Insert the notes before either numpydoc section that can come after 'Notes', otherwise at end
    references_pos = docstring.find("References\n----------")
    examples_pos = docstring.find("Examples\n--------")
    index = min({len(docstring), references_pos, examples_pos} - {-1})

    pieces = [docstring[:index],
              "\n\n" if index == len(docstring) else "",
              "Notes\n-----\n" if existing_notes_pos == -1 else "",
              f"{notes_to_add_data}",
              f"\n\n{docstring[index:]}" if len(docstring) > index else ""]
    return "".join(pieces)


def grid_perimeter(nx, ny):
    """
    Return a sequence of (x, y) grid points for the perimeter of a grid.

    The sequence represents an open path starting at (0, 0); traversing clockwise
    through (0, ny), (nx, ny), and (nx, 0); and then returning to (0, 0)
    exclusive (i.e., stopping at (1, 0)).

    Parameters
    ----------
    nx : `int`
        The number of grid cells in the X direction
    ny : `int`
        The number of grid cells in the Y direction

    Returns
    -------
    `numpy.ndarray`
        A Nx2 array of (x, y) grid points, where N is ``2 * (nx + ny)``
    """
    edges = [[np.zeros(ny), np.arange(ny)],  # left edge
             [np.arange(nx), np.full(nx, ny)],  # top edge
             [np.full(ny, nx), np.arange(ny, 0, -1)],  # right edge
             [np.arange(nx, 0, -1), np.zeros(nx)]]  # bottom edge
    return np.hstack(edges).T


def extent_in_other_wcs(original_wcs, target_wcs, *, method, original_shape=None, integers=False):
    """
    Returns the pixel extent of one WCS in a different WCS.

    Parameters
    ----------
    original_wcs : `~astropy.wcs.WCS`
        The original WCS
    target_wcs : `~astropy.wcs.WCS`
        The target WCS
    method : `str`
        The method for determining the extent: 'all', 'edges', or 'corners'
    original_shape: 2-element tuple
        The array shape of the original WCS.
        This is optional if it is already defined in ``original_wcs``
    integers : `bool`
        If `True`, round the output appropriately to integer values.
        Defaults to `False`.

    Returns
    -------
    min_x, max_x, min_y, max_y : `float` or `int`
        The pixel extent
    """
    ny, nx = original_wcs.array_shape if original_shape is None else original_shape
    if method == 'all':
        pixels = np.indices((nx + 1, ny + 1)).reshape((2, -1)) - 0.5
    elif method == 'edges':
        pixels = grid_perimeter(nx, ny).T - 0.5
    elif method == 'corners':
        pixels = np.array([[0, 0, nx, nx], [0, ny, ny, 0]]) - 0.5
    else:
        raise ValueError("The allowed options for `method` are 'all', 'edges', or 'corners'.")

    xy = np.stack(pixel_to_pixel(original_wcs, target_wcs, *pixels))

    if not np.all(np.isfinite(xy)):
        if method == 'corners':
            raise RuntimeError("The extent could not be automatically determined from the corners. "
                               "Try specifying 'all' or 'edges'.")
        elif method == 'edges':
            raise RuntimeError("The extent could not be automatically determined from the edges. "
                               "Try specifying 'all'.")
        else:
            raise RuntimeError("The extent could not be automatically determined because all of "
                               "the coordinates in the original WCS transformed to NaNs.")

    min_xy = np.min(xy, axis=1)
    max_xy = np.max(xy, axis=1)
    if integers:
        min_xy = (np.floor(min_xy + 0.5)).astype(int)
        max_xy = (np.ceil(max_xy - 0.5)).astype(int)
    return [min_xy[0], max_xy[0], min_xy[1], max_xy[1]]
