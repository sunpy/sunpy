# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
'''
Main class for representing cubes - 3D sets of data where one axis is
a spectral dimension.
'''
# NOTE: This module uses version 1.02 of "Time coordinates in FITS" by
# Rots et al, available at http://hea-www.cfa.harvard.edu/~arots/TimeWCS/
# This draft standard may change.

import numpy as np
import matplotlib.pyplot as plt
import astropy.nddata
from sunpy.map import GenericMap
from sunpy.visualization.imageanimator import ImageAnimator
import astropy.units as u
from astropy.wcs._wcs import InconsistentAxisTypesError
from sunpy.wcs import wcs_util

__all__ = ['Cube', 'CubeError']


class Cube(astropy.nddata.NDData):
    ''' Class representing spectral cubes.

        Attributes
        ----------
        data: numpy ndarray
            The spectral cube holding the actual data in this object. The axes'
            priorities are time, spectral, celestial. This means that if
            present, each of these axis will take precedence over the others.
            For example, in an x, y, t cube the order would be (t,x,y) and in a
            lambda, t, y cube the order will be (t, lambda, y).

        axes_wcs: astropy WCS object
            The WCS object containing the axes' information

        meta: dict
            Header containing the wavelength-specific metadata as well as the
            whole-file metadata
    '''

    def __init__(self, data, wcs, meta=None, **kwargs):
        data, wcs = _orient(data, wcs)
        astropy.nddata.NDData.__init__(self, data=data,
                                       meta=meta,
                                       **kwargs)
        self.axes_wcs = wcs
        # We don't send this to NDData because it's not
        # supported as of astropy 0.3.2. Eventually we will.
        # Also it's called axes_wcs because wcs belongs to astropy.nddata and
        # that messes up slicing.

    def plot_wavelength_slice(self, offset, axes=None,
                              style='imshow', **kwargs):
        '''
        Plots an x-y graph at a certain specified wavelength onto the current
        axes. Keyword arguments are passed on to matplotlib.

        Parameters
        ----------
        offset: int or float
            The offset from the primary wavelength to plot. If it's an int it
            will plot the nth wavelength from the primary; if it's a float then
            it will plot the closest wavelength. If the offset is out of range,
            it will plot the primary wavelength (offset 0)

        axes: matplotlib.axes or None:
            The axes to plot onto. If None the current axes will be used.

        style: 'imshow' or 'pcolormesh'
            The style of plot to be used. Default is 'imshow'
        '''
        if axes is None:
            axes = plt.gca()

        data = self._choose_wavelength_slice(offset)
        if data is None:
            data = self._choose_wavelength_slice(0)

        if style is 'imshow':
            plot = axes.imshow(data, **kwargs)
        elif style is 'pcolormesh':
            plot = axes.pcolormesh(data, **kwargs)

        return plot

    def plot_x_slice(self, offset, axes=None,
                     style='imshow', **kwargs):
        '''
        Plots an x-y graph at a certain specified wavelength onto the current
        axes. Keyword arguments are passed on to matplotlib.

        Parameters
        ----------
        offset: int or float
            The offset from the initial x value to plot. If it's an int it
            will plot slice n from the start; if it's a float then
            it will plot the closest x-distance. If the offset is out of range,
            it will plot the primary wavelength (offset 0)

        axes: matplotlib.axes or None:
            The axes to plot onto. If None the current axes will be used.

        style: 'imshow' or 'pcolormesh'
            The style of plot to be used. Default is 'imshow'
        '''
        if axes is None:
            axes = plt.gca()

        data = self._choose_x_slice(offset)
        if data is None:
            data = self._choose_x_slice(0)

        if style is 'imshow':
            plot = axes.imshow(data, **kwargs)
        elif style is 'pcolormesh':
            plot = axes.pcolormesh(data, **kwargs)

        return plot

    def animate(self, *args, **kwargs):
        '''Plots an interactive visualization of this cube with a slider
        controlling the wavelength axis.
        Parameters other than data are passed to ImageAnimator, which in turn
        passes them to imshow.'''
        i = ImageAnimator(data=self.data, *args, **kwargs)
        return i

    def _choose_wavelength_slice(self, offset):
        '''Retrieves an x-y slice at a wavelength specified by the cube's
        primary wavelength plus the given offset.

        Parameters
        ----------
        offset: int or astropy quantity
            Offset from the cube's primary wavelength. If the value is an int,
            then it returns that slice. Otherwise, it will return the nearest
            wavelength to the one specified.
        '''
        if 'WAVE' not in self.axes_wcs.wcs.ctype:
            raise CubeError(2, "Spectral dimension not present")

        axis = 1 if self.axes_wcs.wcs.ctype[-1] in ['TIME', 'UTC'] else 0
        arr = None
        if (isinstance(offset, int) and offset >= 0 and
           offset < self.data.shape[axis]):
            arr = self.data.take(offset, axis=axis)

        if isinstance(offset, u.Quantity):
            delta = self.axes_wcs.wcs.cdelt[-1 - axis] * u.m
            wloffset = offset.to(u.m) / delta
            wloffset = int(wloffset)
            if wloffset >= 0 and wloffset < self.data.shape[axis]:
                arr = self.data.take(wloffset, axis=axis)

        return arr

    def _choose_x_slice(self, offset):
        '''
        Retrieves a lambda-y slice at an x coordinate specified by the cube's
        primary wavelength plus the given offset.

        Parameters
        ----------
        offset: int or astropy quantity
            Offset from the cube's initial x. If the value is an int,
            then it returns that slice. Otherwise, it will return the nearest
            wavelength to the one specified.
        '''
        arr = None
        axis = 1 if self.axes_wcs.wcs.ctype[-2] != 'WAVE' else 2
        if (isinstance(offset, int) and offset >= 0 and
           offset < self.data.shape[axis]):
            arr = self.data.take(offset, axis=axis)

        if isinstance(offset, u.Quantity):
            unit = self.axes_wcs.wcs.cunit[-1 - axis]
            delta = self.axes_wcs.wcs.cdelt[-1 - axis] * unit
            wloffset = offset.to(unit) / delta
            wloffset = int(wloffset)
            if wloffset >= 0 and wloffset < self.data.shape[axis]:
                arr = self.data.take(wloffset, axis=axis)

        return arr

    def slice_to_map(self, chunk, *args, **kwargs):
        # TODO: implement slice-by-float functionality
        '''
        Converts a given frequency chunk to a SunPy Map. Extra parameters are
        passed on to Map.

        Parameters
        ----------
        chunk: int or float or (int, int) or (float, float)
            The piece of the cube to convert to a map. If it's a single number,
            then it will return that single-slice map, otherwise it will
            aggregate the given range.
        '''
        if self.axes_wcs.wcs.ctype[-2] == 'WAVE':
            raise CubeError(3, "Cannot construct a map with only one spatial dimension")

        if isinstance(chunk, tuple):
            maparray = self.data[chunk[0]:chunk[1], :, :].sum(0)
        else:
            maparray = self.data[chunk, :, :]
        gmap = GenericMap(data=maparray, header=self.meta,
                          *args, **kwargs)
        return gmap

    def slice_to_lightcurve(self, wavelength, y_coord):
        if self.axes_wcs.wcs.ctype[-1] in ['TIME', 'UTC']:
            raise CubeError(1, 'Cannot create a lightcurve with no time axis')
        # TODO: implement this!


def _orient(array, wcs):
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

    order = _select_order(axtypes)
    result_array = array.transpose(order)

    try:
        wcs.get_axis_types()
    except InconsistentAxisTypesError:
        # This means there's an unmatched celestial axis.
        wcs = wcs_util.add_celestial_axis(wcs)
        order = order + [3]

    result_wcs = wcs_util.reindex_wcs(wcs, np.array(order)[::-1])
    return result_array, result_wcs


def _select_order(axtypes):
    order = [(0, t) if t in ['TIME', 'UTC'] else
             (1, t) if t == 'WAVE' else
             (axtypes.index(t) + 2, t) for t in axtypes]
    order.sort()
    result = [axtypes.index(s) for (f, s) in order]
    return result


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
        