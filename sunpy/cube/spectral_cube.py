# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
'''
Main class for representing spectral cubes - 3D sets of data where one axis is
a spectral dimension.
'''

import numpy as np
import matplotlib.pyplot as plt
import astropy.nddata
import sunpy.util.util as util
from sunpy.map import GenericMap
from sunpy.visualization.imageanimator import ImageAnimator
import astropy.units as u

__all__ = ['SpectralCube']


class SpectralCube(astropy.nddata.NDData):
    ''' Class representing spectral cubes.

        Attributes
        ----------
        data: numpy ndarray
            The spectral cube holding the actual data in this object. The axes
            are always [spectral dimension, spatial dimension, extra dimension]
            where the extra dimension can be time or another spatial dimension.

        wcs: astropy WCS object
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
        self.wcs = wcs
        # We don't send this to NDData because it's not
        # supported as of astropy 0.3.2. Eventually we will.

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
            data = self.data[0, :, :]

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
            data = self.data[:, 0, :]

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
        arr = None
        if (isinstance(offset, int) and offset >= 0 and
           offset < len(self.data)):
            arr = self.data[offset, :, :]

        if isinstance(offset, u.Quantity):
            delta = self.wcs.wcs.cdelt[2] * u.m
            wloffset = offset.to(u.m) / delta
            wloffset = int(wloffset)
            if wloffset >= 0 and wloffset < len(self.data):
                arr = self.data[wloffset, :, :]

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
        if (isinstance(offset, int) and offset >= 0 and
           offset < self.data.shape[2]):
            arr = self.data[:, :, offset]
            arr = arr.T

        if isinstance(offset, u.Quantity):
            unit = self.wcs.wcs.cunit[0]
            delta = self.wcs.wcs.cdelt[0] * unit
            wloffset = offset.to(unit) / delta
            wloffset = int(wloffset)
            if wloffset >= 0 and wloffset < self.data.shape[2]:
                arr = self.data[:, :, wloffset]
                arr = arr.T

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
        if isinstance(chunk, tuple):
            maparray = self.data[chunk[0]:chunk[1], :, :].sum(0)
        else:
            maparray = self.data[chunk, :, :]
        gmap = GenericMap(data=np.array(maparray), header=self.meta,
                          *args, **kwargs)
        return gmap


def _orient(array, wcs):
    # This is mostly lifted from astropy's spectral cube.
    """
    Given a 3-d spectral cube and WCS, swap around the axes so that the
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

    # reverse from wcs -> numpy convention
    axtypes = wcs.get_axis_types()[::-1]

    types = [a['coordinate_type'] for a in axtypes]

    # header sanitization removed to allow for arbitrary axes, including time.

    nums = [None if a['coordinate_type'] == 'spectral' else a['number']
            for a in axtypes]

    if 'stokes' in types:
        raise ValueError("Input WCS should not contain stokes")

    order = [types.index('spectral'), nums.index(1), nums.index(0)]
    result_array = array.transpose(order)

    order = wcs.wcs.naxis - np.array(order[::-1]) - 1
    result_wcs = util.reindex_wcs(wcs, order)

    return result_array, result_wcs
    