# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import numpy as np

from matplotlib import pyplot as plt

class Spectrum(np.ndarray):
    """ Class representing a spectrum.
    
    Parameters
    ----------
    data : np.ndarray
        one-dimensional array which the intensity at a particular frequency
        at every data-point.
    freq_axis : np.ndarray
        one-dimensional array with the frequency values at every data point
    """
    def __new__(cls, data, *args, **kwargs):
        return np.asarray(data).view(cls)

    def __init__(self, data, freq_axis):
        self.freq_axis = freq_axis

    def plot(self, overlays=[], **matplotlib_args):
        """
        Plot spectrum onto figure.
        
        Parameters
        ----------
        figure : matplotlib.figure.Figure
            Figure to plot the spectrogram on. If None, new Figure is created.
        overlays : list
            List of overlays (functions that receive figure and axes and return
            new ones) to be applied after drawing.
        """
        # [] as default argument is okay here because it is only read.
        # pylint: disable=W0102,R0914

        if figure is None:
            figure = plt.figure(frameon=True)
        axes = figure.add_subplot(111)

        params = {}

        params.update(matplotlib_args)

        axes.plot(self.freq_axis, self, **params)

        for overlay in overlays:
            figure, axes = overlay(figure, axes)
        return figure