# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
''''
Utilities used in the sunpy.cube.cube module. Moved here to prevent clutter and
aid readability.
'''

import warnings
from copy import deepcopy
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


def handle_slice_to_spectrum(cube, item):
    '''
    Given a cube and a getitem argument, with the knowledge that that slice
    represents a spectrum, return the spectrum that corresponds to that slice.

    Parameters
    ----------
    cube: sunpy.cube.Cube
        The cube to slice
    item: int or slice object or tuple of these
        The slice to make
    '''
    if isinstance(item, int):
        spec = cube.slice_to_spectrum(item, None)
    elif iter_isinstance(item, int, slice, int):
        spec = cube.slice_to_spectrum(item[0], item[2])
    elif iter_isinstance(item, slice, int, int):
        spec = cube.slice_to_spectrum(item[1], item[2])
    else:
        spec = cube.slice_to_spectrum(item[0], None)
    return spec


def handle_slice_to_lightcurve(cube, item):
    '''
    Given a cube and a getitem argument, with the knowledge that that slice
    represents a lightcurve, return the lightcurve that corresponds to that
    slice.

    Parameters
    ----------
    cube: sunpy.cube.Cube
        The cube to slice
    item: int or slice object or tuple of these
        The slice to make
    '''
    if iter_isinstance(item, slice, int, int):
        lightc = cube.slice_to_lightcurve(item[1], item[2])
    else:
        lightc = cube.slice_to_lightcurve(item[1])
    return lightc


def handle_slice_to_map(cube, item):
    '''
    Given a cube and a getitem argument, with the knowledge that that slice
    represents a map, return the map that corresponds to that slice.

    Parameters
    ----------
    cube: sunpy.cube.Cube
        The cube to slice
    item: int or slice object or tuple of these
        The slice to make
    '''
    if isinstance(item, int):
        gmap = cube.slice_to_map(item)
    else:
        gmap = cube.slice_to_map(item[0])
    return gmap


def reduce_dim(cube, axis, keys):
    '''
    Given an axis and a slice object, returns a new cube with the slice
    applied along the given dimension. For example, in a time-x-y cube,
    a reduction along the x axis (axis 1) with a slice value (1, 4, None)
    would return a cube where the only x values were 1 to 3 of the original
    cube.

    Parameters
    ----------
    cube: sunpy.cube.Cube
        The cube to reduce
    axis: int
        The dimension to reduce
    keys: slice object
        The slicing to apply
    '''
    waxis = -1 - axis
    start = keys.start if keys.start is not None else 0
    stop = keys.stop if keys.stop is not None else cube.data.shape[axis]
    if stop > cube.data.shape[axis]:
        stop = cube.data.shape[axis]
    if start < 0:
        start = 0
    step = keys.step if keys.step is not None else 1
    indices = range(start, stop, step)
    newdata = cube.data.take(indices, axis=axis)
    newwcs = cube.axes_wcs.deepcopy()
    if keys.step is not None:
        newwcs.wcs.cdelt[waxis] *= keys.step
    if keys.start is not None:
        start = keys.start
        newwcs.wcs.crpix[waxis] = 0
        newwcs.wcs.crval[waxis] = (cube.axes_wcs.wcs.crval[waxis] +
                                   cube.axes_wcs.wcs.cdelt[waxis] * start)
    newcube = deepcopy(cube)
    newcube.data = newdata
    newcube.axes_wcs = newwcs
    return newcube


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
