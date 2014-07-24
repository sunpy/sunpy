# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
''''
Utilities used in the sunpy.cube.cube module. Moved here to prevent clutter and
aid readability.
'''

import warnings
from astropy.wcs._wcs import InconsistentAxisTypesError
import numpy as np
from sunpy.wcs import wcs_util


def orient(array, wcs):
    # This is mostly lifted from astropy's spectral cube.
    """
    Given a 3-d cube and a WCS, swap around the axes so that the
    spectral axis cube is the first in Numpy notation, and the last in WCS
    notation.

    Parameters
    ----------
    array : `~numpy.ndarray`
        The input 3-d array with two position dimensions and one spectral
        dimension.
    wcs : `~astropy.wcs.WCS`
        The input 3-d WCS with two position dimensions and one spectral
        dimension.
    """

    if array.ndim != 3:
        raise ValueError("Input array must be 3-dimensional")

    if wcs.wcs.naxis != 3:
        raise ValueError("Input WCS must be 3-dimensional")

    axtypes = list(wcs.wcs.ctype)

    array_order = select_order(axtypes[::-1])
    result_array = array.transpose(array_order)

    try:
        wcs.get_axis_types()
    except InconsistentAxisTypesError:
        warnings.warn("Only one spatial axis found. Adding another one...",
                      UserWarning)
        wcs = wcs_util.add_celestial_axis(wcs)

    wcs_order = np.array(select_order(list(wcs.wcs.ctype)))[::-1]
    result_wcs = wcs_util.reindex_wcs(wcs, wcs_order)
    return result_array, result_wcs


def select_order(axtypes):
    '''
    Returns the indices of the correct axis priority for the given list of WCS
    CTYPEs. For example, given ['HPLN-TAN', 'TIME', 'WAVE'] it will return
    [1, 2, 0] because index 1 (time) has the highest priority, followed by
    wavelength and finally solar-x. When two or more celestial axes are in the
    list, order is preserved between them (i.e. only TIME, UTC and WAVE are
    moved)

    Parameters
    ----------
    axtypes: str list
        The list of CTYPEs to be modified.
    '''
    order = [(0, t) if t in ['TIME', 'UTC'] else
             (1, t) if t == 'WAVE' else
             (axtypes.index(t) + 2, t) for t in axtypes]
    order.sort()
    result = [axtypes.index(s) for (_, s) in order]
    return result


def iter_isinstance(obj, *types):
    '''
    Given an iterable object and a list of types, classes or tuples of types
    and classes determine if the given object's items are instances of the
    given types.

    Parameters
    ----------
    obj: tuple
        The object to check
    *types: any number of types or classes
        The classes to check against
    '''
    if not isinstance(obj, tuple) or len(obj) != len(types):
        return False
    return all(isinstance(o, t) for o, t in zip(obj, types))


class CubeError(Exception):
    '''
    Class for handling Cube errors.
    '''
    errors = {0: 'Unspecified error',
              1: 'Time dimension not present',
              2: 'Spectral dimension not present',
              3: 'Insufficient spatial dimensions'}

    def __init__(self, value, msg):
        self.value = value
        self.message = msg

    def __str__(self):
        return 'ERROR ' + repr(self.value) + ' (' \
               + self.errors.get(self.value, '') + '): ' + self.message
