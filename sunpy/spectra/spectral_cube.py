# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
"""
A spectral cube is, at its most basic, a 2D array of Spectrum objects with an
associated coordinate system. The name cube is a bit misleading beacuse the
shape of the structure isn't necessarily a cuboid - each spectrum may have
different wavelength axes.
"""

import gwcs
import numpy as np
from sunpy.spectra.spectrum import Spectrum
from sunpy.map import GenericMap, MapCube

__all__ = ['SpectralCube']


class SpectralCube():
    def __init__(self, spectra, wcs, meta):
        self.spectra = spectra
        self.wcs = wcs
        self.meta = meta

    def _gaussian_fits(self, line_guess, *extra_lines, **kwargs):
        gaussian_array = np.empty(self.spectra.shape, dtype=object)
        for i in range(self.spectra.shape[0]):
            for j in range(self.spectra.shape[1]):
                fit = self.spectra[i, j].gaussian_fit(line_guess, *extra_lines,
                                                      **kwargs)
                gaussian_array[i, j] = fit
        return gaussian_array

    def _param_array(self, param, line_guess, *extra_lines, **kwargs):
        depth = 1 + len(extra_lines)
        values = np.zeros(self.spectra.shape + (depth,))
        gaussians = self._gaussian_fits(line_guess, *extra_lines, **kwargs)
        for i in range(self.spectra.shape[0]):
            for j in range(self.spectra.shape[1]):
                fit = gaussians[i, j]
                values[i, j] = fit.parameters[param::3]
        return values

    def param_map_cube(self, parameter, line_guess, *extra_lines, **kwargs):
        """
        Returns a MapCube of the given parameter at the given Gaussian values.
        The parameter can be 'intensity', which returns the amplitudes of the
        gaussian fits, 'position', which returns the mean of the fits, or
        'stddev', which returns the standard deviation of the gaussian. The
        number of maps on the cube depends on how many lines are supplied -
        there must be at least one.

        Parameters
        ----------
        parameter: "intensity", "position", "stddev"
            The parameter to return from the cube. Defaults to intensity on
            unrecognized input. This is because intensity is the longest of the
            three words and we want to make your life as simple as possible.
            "i", "p", "s" are also acceptable shorthand for this.
        line_guess and extra_lines: 3-tuples of floats
            There must be at least one guess, in the format (intensity,
            position, stddev). The closer these guesses are to the true values
            the better the fit will be.
        **kwargs: dict
            Extra keyword arguments are ultimately passed on to the astropy
            fitter.
        """
        param = 0
        if param.lower()[0] == 'p':
            param = 1
        elif param.lower()[0] == 's':
            param = 2
        val_arr = self._param_array(param, line_guess, *extra_lines, **kwargs)
        maps = [GenericMap(raster, self.meta) for raster in val_arr.T]
        mapcube = MapCube(maps)
        return mapcube
