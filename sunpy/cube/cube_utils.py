# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
"""
Utilities used in the sunpy.cube.cube module. Moved here to prevent clutter and
aid readability.
"""

from copy import deepcopy
import numpy as np
from sunpy.wcs import wcs_util
from astropy import units as u


def orient(array, wcs):
    # This is mostly lifted from astropy's spectral cube.
    """
    Given a 3 or 4D cube and a WCS, swap around the axes so that the
    axes are in the correct order: the first in Numpy notation, and the last
    in WCS notation.

    Parameters
    ----------
    array : `~numpy.ndarray`
        The input 3- or 4-d array with two position dimensions and one spectral
        dimension.
    wcs : `~astropy.wcs.WCS`
        The input 3- or 4-d WCS with two position dimensions and one spectral
        dimension.
    """

    if array.ndim != 3 and array.ndim != 4:
        raise ValueError("Input array must be 3- or 4-dimensional")

    if not ((wcs.wcs.naxis == 3 and array.ndim == 3) or
            (wcs.wcs.naxis == 4 and array.ndim == 4 and not wcs.was_augmented)
            or (wcs.wcs.naxis == 4 and array.ndim == 3 and wcs.was_augmented)):
        raise ValueError("WCS must have the same dimensions as the array")

    axtypes = list(wcs.wcs.ctype)

    if wcs.was_augmented:
        array_order = select_order(axtypes[2::-1])
    else:
        array_order = select_order(axtypes)
    result_array = array.transpose(array_order)

    wcs_order = np.array(select_order(axtypes))[::-1]
    result_wcs = wcs_util.reindex_wcs(wcs, wcs_order)
    return result_array, result_wcs


def select_order(axtypes):
    """
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
    """
    order = [(0, t) if t in ['TIME', 'UTC'] else
             (1, t) if t == 'WAVE' else
             (axtypes.index(t) + 2, t) for t in axtypes]
    order.sort()
    result = [axtypes.index(s) for (_, s) in order]
    return result


def iter_isinstance(obj, *types):
    """
    Given an iterable object and a list of types, classes or tuples of types
    and classes determine if the given object's items are instances of the
    given types.

    Parameters
    ----------
    obj: tuple
        The object to check
    *types: any number of types or classes
        The classes to check against
    """
    if not isinstance(obj, tuple) or len(obj) != len(types):
        return False
    return all(isinstance(o, t) for o, t in zip(obj, types))


def handle_slice_to_spectrum(cube, item):
    """
    Given a cube and a getitem argument, with the knowledge that that slice
    represents a spectrum, return the spectrum that corresponds to that slice.

    Parameters
    ----------
    cube: sunpy.cube.Cube
        The cube to slice
    item: int or slice object or tuple of these
        The slice to make
    """
    if cube.data.ndim == 3:
        if isinstance(item, int):
            spec = cube.slice_to_spectrum(item, None)
        elif iter_isinstance(item, int, slice, int):
            spec = cube.slice_to_spectrum(item[0], item[2])
        elif iter_isinstance(item, slice, int, int):
            spec = cube.slice_to_spectrum(item[1], item[2])
        else:
            spec = cube.slice_to_spectrum(item[0], None)
    else:
        if iter_isinstance(item, int, slice, int, int):
            spec = cube.slice_to_spectrum(item[0], item[2], item[3])
        elif (iter_isinstance(item, int, slice, int, slice) or
              iter_isinstance(item, int, slice, int)):
            spec = cube.slice_to_spectrum(item[0], item[2], None)
        elif iter_isinstance(item, int, slice, slice, int):
            spec = cube.slice_to_spectrum(item[0], None, item[3])
    return spec


def handle_slice_to_lightcurve(cube, item):
    """
    Given a cube and a getitem argument, with the knowledge that that slice
    represents a lightcurve, return the lightcurve that corresponds to that
    slice.

    Parameters
    ----------
    cube: sunpy.cube.Cube
        The cube to slice
    item: int or slice object or tuple of these
        The slice to make
    """
    if cube.data.ndim == 3:
        if iter_isinstance(item, slice, int, int):
            lightc = cube.slice_to_lightcurve(item[1], item[2])
        else:
            lightc = cube.slice_to_lightcurve(item[1])
    else:
        if iter_isinstance(item, slice, int, int, int):
            lightc = cube.slice_to_lightcurve(item[1], item[2], item[3])
        elif iter_isinstance(item, slice, int, slice, int):
            lightc = cube.slice_to_lightcurve(item[1], x_coord=item[3])
        else:
            lightc = cube.slice_to_lightcurve(item[1], y_coord=item[2])
    return lightc


def handle_slice_to_map(cube, item):
    """
    Given a cube and a getitem argument, with the knowledge that that slice
    represents a map, return the map that corresponds to that slice.

    Parameters
    ----------
    cube: sunpy.cube.Cube
        The cube to slice
    item: int or slice object or tuple of these
        The slice to convert
    """
    if cube.data.ndim == 3:
        if isinstance(item, int):
            gmap = cube.slice_to_map(item)
        else:
            gmap = cube.slice_to_map(item[0])
    else:
        gmap = cube.slice_to_map(item[0], item[1])
    return gmap


def handle_slice_to_cube(hypcube, item):
    """
    Given a hypercube and a getitem argument, with the knowledge that the slice
    represents a 3D cube, return the cube that corresponds to that slice.

    Parameters
    ----------
    hypcube: sunpy.cube.Cube
        The 4D hypercube to slice
    item: int or slice, or tuple of these
        The slice to convert
    """
    if isinstance(item, int):
        chunk = item
        axis = 0
    else:
        chunk = [i for i in item if isinstance(i, int)][0]
        axis = item.index(chunk)
    return hypcube.slice_to_cube(axis, chunk)


def reduce_dim(cube, axis, keys):
    """
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
    """
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


def getitem_3d(cube, item):
    """
    Handles Cube's __getitem__ method for 3-dimensional cubes.

    Parameters
    ----------
    cube: sunpy.cube object
        The cube to get the item from
    item: int, slice object, or tuple of these
        The item to get from the cube
    """
    axes = cube.axes_wcs.wcs.ctype
    slice_to_map = (axes[-2] != 'WAVE' and
                    (isinstance(item, int) or
                     iter_isinstance(item, int, slice) or
                     iter_isinstance(item, int, slice, slice)))
    slice_to_spectrum = (((isinstance(item, int) or
                           iter_isinstance(item, int, slice) or
                           iter_isinstance(item, int, slice, slice) or
                           iter_isinstance(item, int, slice, int))
                          and axes[-2] == 'WAVE')
                         or (axes[-1] == 'WAVE' and
                             iter_isinstance(item, slice, int, int)))
    slice_to_spectrogram = (iter_isinstance(item, slice, slice, int)
                            and axes[-2] == 'WAVE')
    slice_to_lightcurve = (axes[-2] == 'WAVE' and
                           (iter_isinstance(item, slice, int, int)
                            or iter_isinstance(item, slice, int)
                            or iter_isinstance(item, slice, int, slice)))
    stay_as_cube = (isinstance(item, slice) or
                    (isinstance(item, tuple) and
                     not any(isinstance(i, int) for i in item)))

    reducedcube = reduce_dim(cube, 0, slice(None, None, None))
    if isinstance(item, tuple):
        for i in range(len(item)):
            if isinstance(item[i], slice):
                reducedcube = reduce_dim(reducedcube, i, item[i])

    if slice_to_map:
        result = handle_slice_to_map(reducedcube, item)
    elif slice_to_spectrum:
        result = handle_slice_to_spectrum(reducedcube, item)
    elif slice_to_spectrogram:
        result = reducedcube.slice_to_spectrogram(item[2])
    elif slice_to_lightcurve:
        result = handle_slice_to_lightcurve(reducedcube, item)
    elif stay_as_cube:
        result = reducedcube
    else:
        result = cube.data[item]
    return result


def getitem_4d(cube, item):
    """
    Handles Cube's __getitem__ method for 4-dimensional hypercubes.

    Parameters
    ----------
    cube: sunpy.cube object
        The cube to get the item from
    item: int, slice object, or tuple of these
        The item to get from the cube
    """
    slice_to_map = (iter_isinstance(item, int, int) or
                    iter_isinstance(item, int, int, slice) or
                    iter_isinstance(item, int, int, slice, slice))
    slice_to_spectrogram = iter_isinstance(item, slice, slice, int, int)
    slice_to_spectrum = (iter_isinstance(item, int, slice, int, int) or
                         iter_isinstance(item, int, slice, int) or
                         iter_isinstance(item, int, slice, int, slice) or
                         iter_isinstance(item, int, slice, slice, int))
    slice_to_cube = (isinstance(item, int) or
                     (isinstance(item, tuple) and
                      len([i for i in item if isinstance(i, int)]) == 1))
    slice_to_lightcurve = (iter_isinstance(item, slice, int, int, int) or
                           iter_isinstance(item, slice, int, int) or
                           iter_isinstance(item, slice, int, int, slice) or
                           iter_isinstance(item, slice, int, slice, int))
    stay_as_hypercube = (isinstance(item, slice) or
                         (isinstance(item, tuple) and
                          not any(isinstance(i, int) for i in item)))
    reducedcube = reduce_dim(cube, 0, slice(None, None, None))
    if isinstance(item, tuple):
        for i in range(len(item)):
            if isinstance(item[i], slice):
                reducedcube = reduce_dim(reducedcube, i, item[i])

    if slice_to_map:
        result = handle_slice_to_map(reducedcube, item)
    elif slice_to_spectrum:
        result = handle_slice_to_spectrum(reducedcube, item)
    elif slice_to_spectrogram:
        result = reducedcube.slice_to_spectrogram(item[2], item[3])
    elif slice_to_lightcurve:
        result = handle_slice_to_lightcurve(reducedcube, item)
    elif slice_to_cube:
        result = handle_slice_to_cube(reducedcube, item)
    elif stay_as_hypercube:
        result = reducedcube
    else:
        result = cube.data[item]
    return result


def pixelize_slice(item, wcs):
    """
    Given a getitem slice that may or may not contain astropy units and a wcs,
    convert these to pixels. Raises a CubeError if the units don't match.
    This assumes that the array is not rotated.

    Parameters
    ----------
    item: int, astropy.Quantity, slice, or tuple of these
        The slice to convert to pixels.
    wcs: Sunpy.wcs.wcs.WCS
        The WCS object representing the array/
    """
    if isinstance(item, tuple):
        result = list(range(len(item)))
        for axis in range(len(item)):
            if isinstance(item[axis], slice):
                result[axis] = _convert_slice(item[axis], wcs, axis)
            elif isinstance(item[axis], u.Quantity):
                result[axis] = _convert_point(item[axis].value,
                                              item[axis].unit, wcs, axis)
            else:
                result[axis] = item[axis]
        result = tuple(result)
    elif isinstance(item, u.Quantity):
        result = _convert_point(item.value, item.unit, wcs, 0)
    else:
        result = item

    return result


def _convert_point(value, unit, wcs, axis):
    naxis = wcs.wcs.naxis if not wcs.was_augmented else wcs.wcs.naxis - 1
    cunit = u.Unit(wcs.wcs.cunit[naxis - 1 - axis])
    crpix = wcs.wcs.crpix[naxis - 1 - axis]
    crval = wcs.wcs.crval[naxis - 1 - axis] * cunit
    cdelt = wcs.wcs.cdelt[naxis - 1 - axis] * cunit

    point = (value * unit).to(cunit)
    pointdelta = ((point - crval) / cdelt).value
    point = crpix + pointdelta

    return int(np.round(point))


def _convert_slice(sl, wcs, axis):
    naxis = wcs.wcs.naxis if not wcs.was_augmented else wcs.wcs.naxis - 1
    steps = [sl.start, sl.stop, sl.step]
    values = [None, None, None]
    unit = None
    for i in range(3):
        if isinstance(steps[i], u.Quantity):
            if unit is not None and steps[i].unit != unit:
                raise CubeError(5, "Only one unit per axis may be given")
            else:
                unit = steps[i].unit
                values[i] = steps[i].value
        else:
            values[i] = steps[i]
    if unit is None:
        return sl

    if values[2] is None:
        delta = None
    else:
        cunit = u.Unit(wcs.wcs.cunit[naxis - 1 - axis])
        cdelt = wcs.wcs.cdelt[naxis - 1 - axis] * cunit
        delta = int(np.round(((values[2] * unit).to(cunit) / cdelt).value))

    if values[0] is None:
        start = None
    else:
        start = _convert_point(values[0], unit, wcs, axis)

    if values[1] is None:
        end = None
    else:
        end = _convert_point(values[1], unit, wcs, axis)

    return slice(start, end, delta)


class CubeError(Exception):
    """
    Class for handling Cube errors.
    """
    errors = {0: 'Unspecified error',
              1: 'Time dimension not present',
              2: 'Spectral dimension not present',
              3: 'Insufficient spatial dimensions',
              4: 'Dimension error',
              5: 'Slicing unit error'}

    def __init__(self, value, msg):
        Exception.__init__(self, msg)
        self.value = value

    def __str__(self):
        return 'ERROR ' + repr(self.value) + ' (' \
               + self.errors.get(self.value, '') + '): ' + self.message
