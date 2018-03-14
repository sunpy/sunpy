# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import numpy as np
from matplotlib import pyplot as plt

from sunpy.util.decorators import deprecated

__all__ = ['Spectrum']


@deprecated("0.8", "this module is being moved to radiospectra - http://docs.sunpy.org/projects/radiospectra/")
class Spectrum(np.ndarray):
    """
    Class representing a 1 dimensional spectrum.

    Attributes
    ----------
    freq_axis : `~numpy.ndarray`
        one-dimensional array with the frequency values.

    data\ : `numpy.ndarray`
        One-dimensional array which the intensity at a particular frequency at
        every data-point.


    Examples
    --------
    >>> from sunpy.spectra.spectrum import Spectrum
    >>> import numpy as np
    >>> data = np.linspace(1, 100, 100)
    >>> freq_axis = np.linspace(0, 10, 100)
    >>> spec = Spectrum(data, freq_axis)
    >>> spec.peek()   # doctest: +SKIP
    """
    def __new__(cls, data, *args, **kwargs):
        return np.asarray(data).view(cls)

    def __init__(self, data, freq_axis):
        if np.shape(data)[0] != np.shape(freq_axis)[0]:
            raise ValueError('Dimensions of data and frequency axis do not match')
        self.freq_axis = freq_axis

    def plot(self, axes=None, **matplot_args):
        """
        Plot spectrum onto current axes.

        Parameters
        ----------
        axes : `~matplotlib.axes.Axes` or None
            If provided the spectrum will be plotted on the given axes.
            Else the current matplotlib axes will be used.

        **matplot_args : dict
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        newaxes : `~matplotlib.axes.Axes`
            The plot axes.
        """

        # Get current axes
        if not axes:
            axes = plt.gca()

        params = {}
        params.update(matplot_args)

        lines = axes.plot(self.freq_axis, self, **params)

        return lines

    def peek(self, **matplot_args):
        """
        Plot spectrum onto a new figure. An example is shown below.

        .. plot::

            from sunpy.spectra.spectrum import Spectrum
            import numpy as np
            spec = Spectrum(np.linspace(1, 100, 100), np.linspace(0, 10, 100))
            spec.peek()

        Parameters
        ----------
        **matplot_args : dict
            Any additional plot arguments that should be used
            when plotting.
        """

        figure = plt.figure()
        lines = self.plot(**matplot_args)
        figure.show()
