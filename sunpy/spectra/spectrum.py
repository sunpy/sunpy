# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import numpy as np

from matplotlib import pyplot as plt

__all__ = ['Spectrum']

class Spectrum(np.ndarray):
    """
    Class representing a spectrum.
    
    Attributes
    ----------
    freq_axis : np.ndarray
        one-dimensional array with the frequency values at every data point

    data\ : np.ndarray
        one-dimensional array which the intensity at a particular frequency at every data-point.
    """
    def __new__(cls, data, *args, **kwargs):
        return np.asarray(data).view(cls)

    def __init__(self, data, freq_axis):
        self.data = data
        self.freq_axis = freq_axis

    def plot(self, axes=None, **matplot_args):
        """
        Plot spectrum onto current axes. Behaves like matplotlib.pylot.plot()
        
        Parameters
        ----------
        axes: matplotlib.axes object or None
            If provided the spectrum will be plotted on the given axes. 
            Else the current matplotlib axes will be used.
        """
        
        #Get current axes
        if not axes:
            axes = plt.gca()
        
        params = {}
        params.update(matplot_args)
        
        #This is taken from mpl.pyplot.plot() as we are trying to
        #replicate that functionality
        
        # allow callers to override the hold state by passing hold=True|False
        washold = axes.ishold()
        hold = matplot_args.pop('hold', None)
        
        if hold is not None:
            axes.hold(hold)
        try:
            lines = axes.plot(self.freq_axis, self, **params)
        finally:
            axes.hold(washold)

        return lines
    
    def peek(self, **matplot_args):
        """
        Plot spectrum onto a new figure.
        """
        
        figure = plt.figure()
        
        lines = self.plot(**matplot_args)
        
        figure.show()
        
        return figure

