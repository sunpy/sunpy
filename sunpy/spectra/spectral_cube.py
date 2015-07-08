# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
"""
A spectral cube is, at its most basic, a 2D array of Spectrum objects with an
associated coordinate system. The name cube is a bit misleading beacuse the
shape of the structure isn't necessarily a cuboid - each spectrum may have
different wavelength axes.
"""

import gwcs
import numpy as np
from sunpy.map import GenericMap, MapCube

__all__ = ['SpectralCube']


class SpectralCube(object):
    """
    Class defining spectral cubes: 2-dimensional arrays of Spectrum objects.
    The individual spectra may have different axis so the shape may not
    necessarily be a perfect cuboid.

    Attributes
    ----------
    spectra: Numpy ndarray of Spectrum objects
        The main data held by the cube. There are two axes, and their
        priorities are the same as in Cubes: the dimensions are (time, solar_y)
        or (solar_x, solar_y), depending on the underlying data
    wcs: sunpy.wcs.WCS object
        WCS system describing the two-axis system. The information about the
        spectral axis is stored in the individual Spectrum objects.
    meta: dict
        Metadata for the current mission and observation.
    """

    def __init__(self, spectra, wcs, meta):
        self.spectra = spectra
        self.wcs = wcs
        self.meta = meta
        self._memo = {}

    def _gaussian_fits(self, line_guess, *extra_lines, **kwargs):
        """
        Returns an array of fit objects from which parameters can be extracted
        corresponding to the line guesses provided.

        Parameters
        ----------
        line_guess and extra_lines: 3-tuples of floats
            There must be at least one guess, in the format (intensity,
            position, stddev). The closer these guesses are to the true values
            the better the fit will be.
        recalc=False: boolean
            If True, the gaussian fits will be recalculated, even if there's an
            existing fit for the given wavelengths already in the memo. This
            keyword should be set to True if changing the amplitude or width of
            the fit.
        **kwargs: dict
            Extra keyword arguments are ultimately passed on to the astropy
            fitter.
        """
        recalc = kwargs.get('recalculate', False)
        key = [guess[1] for guess in (line_guess,) + extra_lines]
        if recalc or key not in self._memo:
            gaussian_array = np.empty(self.spectra.shape, dtype=object)
            for i in range(self.spectra.shape[0]):
                for j in range(self.spectra.shape[1]):
                    fit = self.spectra[i, j].gaussian_fit(line_guess,
                                                          *extra_lines,
                                                          **kwargs)
                    gaussian_array[i, j] = fit
            self._memo[key] = gaussian_array
            return gaussian_array
        else:
            return self._memo[key]

    def _param_array(self, param, line_guess, *extra_lines, **kwargs):
        """
        Returns the values of the parameter specified for the fit array. The
        parameter can be 0 for intensity, 1 for mean or 2 for width. Other
        values may throw exceptions or return meaningless values.

        Parameters
        ----------
        param: int, one of (0, 1, 2)
            The parameter to return from the cube. Intensity is 0, mean is 1
            and standard deviation is 2.
        line_guess and extra_lines: 3-tuples of floats
            There must be at least one guess, in the format (intensity,
            position, stddev). The closer these guesses are to the true values
            the better the fit will be.
        recalc=False: boolean
            If True, the gaussian fits will be recalculated, even if there's an
            existing fit for the given wavelengths already in the memo. This
            keyword should be set to True if changing the amplitude or width of
            the fit.
        **kwargs: dict
            Extra keyword arguments are ultimately passed on to the astropy
            fitter.
        """
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
        recalc=False: boolean
            If True, the gaussian fits will be recalculated, even if there's an
            existing fit for the given wavelengths already in the memo. This
            keyword should be set to True if changing the amplitude or width of
            the fit.
        **kwargs: dict
            Extra keyword arguments are ultimately passed on to the astropy
            fitter.
        """
        param = 0
        if parameter.lower()[0] == 'p':
            param = 1
        elif parameter.lower()[0] == 's':
            param = 2
        val_arr = self._param_array(param, line_guess, *extra_lines, **kwargs)
        maps = [GenericMap(raster, self.meta) for raster in val_arr.T]
        mapcube = MapCube(maps)
        return mapcube

    # TODO: __getitem__
