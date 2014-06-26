# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>

import spectral_cube as sc
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.nddata
from sunpy.map import GenericMap

__all__ = ['SpectralCube']


class SpectralCube(astropy.nddata.NDData):
    # TODO: Should the headers be header objects or dicts? I think dicts...
    ''' Class representing spectral cubes.

        Attributes
        ----------
        cube: spectral_cube
            the spectral cube holding data and coordinates

        data_header: astropy.io.fits.Header object
            Header containing the wavelength-specific metadata

        primary_header: astropy.io.fits.Header object
            Main header containing metadata pertaining to the whole file
    '''

    def __init__(self, cube, data_header=None, primary_header=None, **kwargs):
        astropy.nddata.NDData.__init__(data=cube.unmasked_data,
                                       wcs=cube.wcs,
                                       **kwargs)
        self.cube = cube
        self.data_header = data_header
        self.primary_header = primary_header

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
            data = self.cube.unmasked_data[0, :, :]

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
            data = self.cube.unmasked_data[:, 0, :]

        if style is 'imshow':
            plot = axes.imshow(data, **kwargs)
        elif style is 'pcolormesh':
            plot = axes.pcolormesh(data, **kwargs)

        return plot

    def _choose_wavelength_slice(self, offset):
        '''Retrieves an x-y slice at a wavelength specified by the cube's
        primary wavelength plus the given offset.

        Parameters
        ----------
        offset: int or float
            Offset from the cube's primary wavelength. If the value is an int,
            then it returns that slice. Otherwise, it will return the nearest
            wavelength to the one specified.
        '''
        a = None
        if (isinstance(offset, int) and offset >= 0 and
            offset < len(self.cube.spectral_axis)):
            a = self.cube.unmasked_data[offset, :, :]

        # TODO: this currently fails because delta is a numpy vector
        if isinstance(offset, float):
            delta = self.cube.spectral_axis[1] - self.cube.spectral_axis[0]
            wloffset = offset / float(delta.decompose())
            wloffset = int(wloffset)
            if wloffset >= 0 and wloffset < len(self.cube.spectral_axis):
                a = self.cube.unmasked_data[wloffset, :, :]

        return np.array(a)

    def _choose_x_slice(self, offset):
        '''
        Retrieves a lambda-y slice at an x coordinate specified by the cube's
        primary wavelength plus the given offset.

        Parameters
        ----------
        offset: int or float
            Offset from the cube's initial x. If the value is an int,
            then it returns that slice. Otherwise, it will return the nearest
            wavelength to the one specified.
        '''
        a = None
        if (isinstance(offset, int) and offset >= 0 and
            offset < self.cube.shape[2]):
            a = self.cube.unmasked_data[:, :, offset]

        # TODO: This fails because delta is not a scalar (and it actually gets
        # a wavelength slice, but nevermind...)
        if isinstance(offset, float):
            delta = self.cube.spectral_axis[1] - self.cube.spectral_axis[0]
            wloffset = offset / delta
            wloffset = int(wloffset)
            if wloffset >= 0 and wloffset < len(self.cube.spectral_axis):
                a = self.cube.unmasked_data[wloffset, :, :]

        return np.array(a).T

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
            maparray = self.cube.unmasked_data[chunk[0]:chunk[1], :, :].sum(0)
        else:
            maparray = self.cube.unmasked_data[chunk, :, :]
        # TODO: send a proper header
        m = GenericMap(data=maparray, header=None)
        return m
