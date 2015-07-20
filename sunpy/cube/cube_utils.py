# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101, C0330
"""
Utilities used in the sunpy.cube.cube module. Moved here to prevent clutter and
aid readability.
"""

import numpy as np
from sunpy.wcs import wcs_util
from astropy import units as u
import sunpy.cube.cube


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

    if wcs.oriented:  # If this wcs has already been oriented.
        return array, wcs

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
    result_wcs.was_augmented = wcs.was_augmented
    result_wcs.oriented = True
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


def iter_isinstance(obj, *type_tuples):
    """
    Given an iterable object and a list of tuples of types, classes or tuples
    of types and classes determine if the given object's items are instances of
    the given types. iter_isinstance(obj, types_1, types_2) is shorthand
    for iter_isinstance(obj, types_1) or iter_isinstance(obj, types_2).

    Parameters
    ----------
    obj: tuple
        The object to check
    *types: any number of types or classes
        The classes to check against
    """
    result = False
    if not isinstance(obj, (tuple, list)):
        return False
    for types in type_tuples:
        if len(obj) != len(types):
            continue
        result |= all(isinstance(o, t) for o, t in zip(obj, types))
    return result


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
        elif iter_isinstance(item, (int, slice, int)):
            spec = cube.slice_to_spectrum(item[0], item[2])
        elif iter_isinstance(item, (slice, int, int)):
            spec = cube.slice_to_spectrum(item[1], item[2])
        else:
            spec = cube.slice_to_spectrum(item[0], None)
    else:
        if iter_isinstance(item, (int, slice, int, int)):
            spec = cube.slice_to_spectrum(item[0], item[2], item[3])
        elif iter_isinstance(item, (int, slice, int, slice),
                                   (int, slice, int)):
            spec = cube.slice_to_spectrum(item[0], item[2], None)
        elif iter_isinstance(item, (int, slice, slice, int)):
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
        if iter_isinstance(item, (slice, int, int)):
            lightc = cube.slice_to_lightcurve(item[1], item[2])
        else:
            lightc = cube.slice_to_lightcurve(item[1])
    else:
        if iter_isinstance(item, (slice, int, int, int)):
            lightc = cube.slice_to_lightcurve(item[1], item[2], item[3])
        elif iter_isinstance(item, (slice, int, slice, int)):
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
    meta = cube.meta
    unit = cube.unit
    uncertainty = cube.uncertainty
    mask = cube.mask

    newcube = sunpy.cube.Cube(data=newdata, wcs=newwcs, meta=meta, unit=unit,
                              uncertainty=uncertainty, mask=mask)
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
                     iter_isinstance(item, (int, slice), (int, slice, slice))))
    slice_to_spectrum = (((isinstance(item, int) or
                           iter_isinstance(item, (int, slice),
                                                 (int, slice, slice),
                                                 (int, slice, int)))
                          and axes[-2] == 'WAVE')
                         or (axes[-1] == 'WAVE' and
                             iter_isinstance(item, (slice, int, int))))
    slice_to_spectrogram = (iter_isinstance(item, (slice, slice, int)) and
                            axes[-2] == 'WAVE')
    slice_to_lightcurve = (axes[-2] == 'WAVE' and
                           (iter_isinstance(item, (slice, int, int),
                                                  (slice, int),
                                                  (slice, int, slice))))
    stay_as_cube = (isinstance(item, slice) or
                    (isinstance(item, tuple) and
                     not any(isinstance(i, int) for i in item)))

    reducedcube = reduce_dim(cube, 0, slice(None, None, None))
    # We're not actually reducing a cube, just a way of copying the cube.
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
    slice_to_map = iter_isinstance(item, (int, int), (int, int, slice),
                                         (int, int, slice, slice))
    slice_to_spectrogram = iter_isinstance(item, (slice, slice, int, int))
    slice_to_spectrum = iter_isinstance(item, (int, slice, int, int),
                                              (int, slice, int),
                                              (int, slice, int, slice),
                                              (int, slice, slice, int))
    slice_to_cube = (isinstance(item, int) or
                     (isinstance(item, tuple) and
                      len([i for i in item if isinstance(i, int)]) == 1))
    slice_to_lightcurve = iter_isinstance(item, (slice, int, int, int),
                                                (slice, int, int),
                                                (slice, int, int, slice),
                                                (slice, int, slice, int))
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


def pixelize_slice(item, wcs, _source='cube'):
    """
    Given a getitem slice that may or may not contain astropy units and a wcs,
    convert these to pixels. Raises a CubeError if the units don't match.
    This assumes that the array is not rotated.

    Parameters
    ----------
    item: int, astropy.Quantity, slice, or tuple of these
        The slice to convert to pixels.
    wcs: Sunpy.wcs.wcs.WCS
        The WCS object representing the array
    """
    if isinstance(item, tuple):
        result = list(range(len(item)))
        for axis in range(len(item)):
            if isinstance(item[axis], slice):
                result[axis] = _convert_slice(item[axis], wcs,
                                              axis, _source=_source)
            elif isinstance(item[axis], u.Quantity):
                result[axis] = convert_point(item[axis].value,
                                             item[axis].unit, wcs, axis,
                                             _source=_source)
            else:
                result[axis] = item[axis]
        result = tuple(result)
    elif isinstance(item, u.Quantity):
        result = convert_point(item.value, item.unit, wcs, 0)
    else:
        result = item

    return result


def convert_point(value, unit, wcs, axis, _source='cube'):
    """
    Takes a point on an axis specified by the given wcs and returns the pixel
    coordinate.

    Parameters
    ----------
    value: int or float
        The magnitude of the specified point.
    unit: astropy.unit.Unit
        The unit for the given value. Note this doesn't take in a quantity to
        simplify _convert_slice.
    wcs: sunpy.wcs.wcs.WCS
        The WCS describing the axes system.
    axis: int
        The axis the value corresponds to, in numpy-style ordering (i.e.
        opposite WCS convention)
    """
    if value is None:
        return None  # This is used to simplify None coordinates during slicing
    if unit is None or unit == u.pix or unit == u.pixel:
        return int(value)
    if isinstance(value, u.Quantity):
        value = value.value
        unit = value.unit
    if _source is 'cube':
        wcsaxis = -1 - axis if wcs.oriented or not wcs.was_augmented \
                  else -2 - axis
    else:
        wcsaxis = 1 - axis
    cunit = u.Unit(wcs.wcs.cunit[wcsaxis])
    crpix = wcs.wcs.crpix[wcsaxis]
    crval = wcs.wcs.crval[wcsaxis] * cunit
    cdelt = wcs.wcs.cdelt[wcsaxis] * cunit

    point = (value * unit).to(cunit)
    pointdelta = ((point - crval) / cdelt).value
    point = crpix + pointdelta

    return int(np.round(point))


def _convert_slice(item, wcs, axis, _source='cube'):
    """
    Takes in a slice object that may or may not contain units and translates it
    to pixel coordinates along the given wcs and axis. If there are no units,
    returns the same slice; if there is more than one non-identical unit,
    raises a CubeError.

    Parameters
    ----------
    item: slice object
        The slice to convert to pixels. It may be composed of Nones, ints,
        floats or astropy Quantities, in any combination (with the restriction
        noted above)
    wcs: sunpy.wcs.wcs.WCS
        The WCS describing this system
    axis: int
        The axis the slice corresponds to, in numpy-style ordering (i.e.
        opposite WCS convention)
    """
    if _source is 'cube':
        wcs_ax = -2 - axis if wcs.was_augmented and not wcs.oriented \
                else -1 - axis
    else:
        wcs_ax = 1 - axis
    steps = [item.start, item.stop, item.step]
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
        return item

    if values[2] is None:
        delta = None
    else:
        cunit = u.Unit(wcs.wcs.cunit[wcs_ax])
        cdelt = wcs.wcs.cdelt[wcs_ax] * cunit
        delta = int(np.round(((values[2] * unit).to(cunit) / cdelt).value))

    if values[0] is None:
        start = None
    else:
        start = convert_point(values[0], unit, wcs, axis, _source=_source)

    if values[1] is None:
        end = None
    else:
        end = convert_point(values[1], unit, wcs, axis, _source=_source)

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
              5: 'Slicing unit error',
              6: 'Unaligned array error'}

    def __init__(self, value, msg):
        Exception.__init__(self, msg)
        self.value = value

    def __str__(self):
        return 'ERROR ' + repr(self.value) + ' (' \
               + self.errors.get(self.value, '') + '): ' + self.message
