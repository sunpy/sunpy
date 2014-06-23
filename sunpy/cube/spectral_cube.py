# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>

import spectral_cube as sc
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.nddata


class SpectralCube(astropy.nddata.NDData):
    ''' Class representing spectral cubes.

        Attributes
        ----------
        cube: spectral_cube
            the spectral cube holding data and coordinates

        data_header: astropy.io.fits.Header
            Header containing the wavelength-specific metadata

        primary_header: astropy.io.fits.Header
            Main header containing metadata pertaining to the whole file
    '''

    def __init__(self, cube, data_header=None, primary_header=None):
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

        data = SpectralCube._choose_wavelength_slice(offset)
        if data is None:
            data = self.cube.unmasked_data[0, :, :]

        if style is 'imshow':
            plot = axes.imshow(data, **kwargs)
        elif style is 'pcolormesh':
            plot = axes.pcolormesh(data, **kwargs)

        return plot

    def _choose_cube(self, wavelength):
        '''Retrieves the cube at the specified wavelength.

        Parameters
        ----------
        wavelength: str or int
            If str, the wavelength or spectral line to be fetched. If int, the
            index to be returned.
        '''
        if ((isinstance(wavelength, int) and
             wavelength >= 0 and wavelength < len(self.data_cubes))):
            return self.data_cubes.values()[wavelength]

        if (isinstance(wavelength, str)):
            return self.data_cubes.get(wavelength)

        return None

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
        if (isinstance(offset, int) and offset >= 0 and
            offset < len(self.cube.spectral_axis)):
            return self.cube.unmasked_data[offset, :, :]

        if isinstance(offset, float):
            delta = self.cube.spectral_axis[1] - self.cube.spectral_axis[0]
            wloffset = offset / delta
            wloffset = int(wloffset)
            if wloffset >= 0 and wloffset < len(self.cube.spectral_axis):
                return self.cube.unmasked_data[wloffset, :, :]

            return None
