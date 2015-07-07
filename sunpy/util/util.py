from __future__ import absolute_import

import os
import types
import warnings
from itertools import izip, imap, count

import numpy as np

__all__ = ['to_signed', 'unique', 'print_table',
           'replacement_filename', 'goes_flare_class', 'merge', 'common_base',
           'minimal_pairs', 'polyfun_at', 'expand_list',
           'expand_list_generator', 'Deprecated', 'savitzky_golay']


def to_signed(dtype):
    """ Return dtype that can hold data of passed dtype but is signed.
    Raise ValueError if no such dtype exists.

    Parameters
    ----------
    dtype : np.dtype
        dtype whose values the new dtype needs to be able to represent.
    """
    if dtype.kind == "u":
        if dtype.itemsize == 8:
            raise ValueError("Cannot losslessy convert uint64 to int.")
        dtype = "int{0:d}".format(min(dtype.itemsize * 2 * 8, 64))
    return np.dtype(dtype)


def goes_flare_class(gcls):
    """Convert GOES classes into a number to aid size comparison.  Units are
    watts per meter squared."""
    def calc(gcls):
        powers_of_ten = {'A': 1e-08, 'B': 1e-07, 'C': 1e-06,
                         'M': 1e-05, 'X': 1e-04}
        power = gcls[0].upper()
        if power in powers_of_ten:
            return powers_of_ten[power] * float(gcls[1:])
        else:
            return None

    if isinstance(gcls, types.StringType):
        return calc(gcls)
    if isinstance(gcls, types.ListType):
        return [calc(x) for x in gcls]


def unique(itr, key=None):
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
    width = [max(imap(len, col)) for col in izip(*lst)]
    return linesep.join(
        colsep.join(
            col.ljust(n) for n, col in izip(width, row)
        ) for row in lst
    )


def findpeaks(a):
    """ Find local maxima in 1D. Use findpeaks(-a) for minima. """
    return np.nonzero((a[1:-1] > a[:-2]) & (a[1:-1] > a[2:]))[0]


def polyfun_at(coeff, p):
    """ Return value of polynomial with coefficients (highest first) at
    point (can also be an np.ndarray for more than one point) p. """
    return np.sum(k * p ** n for n, k in enumerate(reversed(coeff)))


def minimal_pairs(one, other):
    """ Find pairs of values in one and other with minimal distance.
    Assumes one and other are sorted in the same sort sequence.

    one, other : sequence
        Sequence of scalars to find pairs from.
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
    """ Given two sorted sequences one and other, for every element
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
    """ Find class that every item of objs is an instance of. """
    for cls in objs[0].__class__.__mro__:
        if all(isinstance(obj, cls) for obj in objs):
            break
    return cls


def merge(items, key=(lambda x: x)):
    """ Given sorted lists of iterables, return new iterable that returns
    elements of all iterables sorted with respect to key. """
    state = {}
    for item in map(iter, items):
        try:
            first = item.next()
        except StopIteration:
            continue
        else:
            state[item] = (first, key(first))

    while state:
        for item, (value, tk) in state.iteritems():
            # Value is biggest.
            if (all(tk >= k for it, (v, k)
               in state.iteritems() if it is not item)):
                yield value
                break
        try:
            n = item.next()
            state[item] = (n, key(n))
        except StopIteration:
            del state[item]


def replacement_filename(path):
    """ Return replacement path for already used path. Enumerates
    until an unused filename is found. E.g., "/home/florian/foo.fits"
    becomes "/home/florian/foo.0.fits", if that is used
    "/home/florian/foo.1.fits", etc. """
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


# =============================================================================
# expand list from :http://stackoverflow.com/a/2185971/2486799
# =============================================================================
def expand_list(input):
    return [item for item in expand_list_generator(input)]


def expand_list_generator(input):
    for item in input:
        if type(item) in [list, tuple]:
            for nested_item in expand_list_generator(item):
                yield nested_item
        else:
            yield item


# =============================================================================
# Deprecation decorator: http://code.activestate.com/recipes/391367-deprecated/
# and http://www.artima.com/weblogs/viewpost.jsp?thread=240845
# =============================================================================
class Deprecated(object):
    """ Use this decorator to deprecate a function or method, you can pass an
    additional message to the decorator:

    @Deprecated("no more")
    """
    def __init__(self, message=""):
        self.message = message

    def __call__(self, func):
        def newFunc(*args, **kwargs):
            warnings.warn("Call to deprecated function {0}. \n {1}".format(
                                                                func.__name__,
                                                                self.message),
                          category=Warning, stacklevel=2)
            return func(*args, **kwargs)

        newFunc.__name__ = func.__name__
        newFunc.__doc__ = func.__doc__
        newFunc.__dict__.update(func.__dict__)
        return newFunc


# =============================================================================
# Savitzky-Golay Filter for smoothing from:
# http://wiki.scipy.org/Cookbook/SavitzkyGolay
# =============================================================================
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.

    Parameters
    ----------
    y : array_like, shape (N,)
        The values of the time history of the signal.
    window_size : int
        The length of the window. Must be an odd integer number.
    order : int
        The order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        The order of the derivative to compute.
        Default of 0 means only smoothing
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.

    Examples
    --------
    >>> t = np.linspace(-4, 4, 500)
    >>> y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    >>> ysg = savitzky_golay(y, window_size=31, order=4)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(t, y, label='Noisy signal')
    >>> plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    >>> plt.plot(t, ysg, 'r', label='Filtered signal')
    >>> plt.legend()
    >>> plt.show()

    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window,
                                                           half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')
