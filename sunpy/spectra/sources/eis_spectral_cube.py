# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>

from sunpy.spectra.spectral_cube import SpectralCube
from sunpy.cube import cube_utils as cu
from sunpy.util.util import savitzky_golay
import numpy as np


class EISSpectralCube(SpectralCube):
    # TODO: init, docstrings, etc.
    def orbital_correction(self, line_guess, ymin=None, ymax=None, **kwargs):
        if not ymin:
            ymin = 0
        else:
            ymin = cu.convert_point(ymin.value, ymin.unit, self.wcs, 1)
        if not ymax:
            ymax = self.spectra.shape[1]
        else:
            ymax = cu.convert_point(ymax.value, ymax.unit, self.wcs, 1)
        if ymin > ymax:
            temp = ymin
            ymin = ymax
            ymax = temp
        # First, get the centers of the line we're correcting for, throughout
        # the whole cube. This is a 2D array.
        centers = self._param_array(1, line_guess, **kwargs)
        # Now get a 1-D array of the average intensities along the y-axis.
        # Ideally the range should be chosen so that this is quiet sun.
        averages = [np.average(arr[ymin:ymax]) for arr in centers]
        # Calculate the offset from the centroid of the emission line over the
        # given data range.
        corrections = line_guess[1] - averages
        # Remove noise by appplying a smoothing filter and getting rid of
        # the less important frequencies.
        window_size = int(corrections.shape[0] / 10)
        window_size += 1 if window_size % 2 == 0 else 0
        smooth_averages = savitzky_golay(corrections, window_size, 3)
        # Apply corrections to the spectra.
        for arr, corr in zip(self.spectra.T, smooth_averages):
            for spec in arr:
                spec.shift_axis(corr)
