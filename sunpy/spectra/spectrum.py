# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>,
#         Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101

from __future__ import absolute_import

import astropy.nddata
from astropy.modeling import models, fitting

from matplotlib import pyplot as plt

__all__ = ['Spectrum']


class Spectrum(astropy.nddata.NDDataArray):
    """
    Class representing a spectrum.

    Attributes
    ----------
    freq_axis : np.ndarray
        one-dimensional array with the frequency values at every data point

    data : np.ndarray
        one-dimensional array which the intensity at a particular frequency at
        every data-point.
    """

    def __init__(self, data, freq_axis, **kwargs):
        astropy.nddata.NDDataArray.__init__(self, data=data, **kwargs)
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

        # Get current axes
        if not axes:
            axes = plt.gca()

        params = {}
        params.update(matplot_args)

        # This is taken from mpl.pyplot.plot() as we are trying to
        # replicate that functionality

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
        self.plot(**matplot_args)
        figure.show()
        return figure

    def shift_axis(self, offset):
        """
        Shifts the entire wavelength axis by a given linear offset

        Parameters
        ----------
        offset: float
            The amount to offset by
        """
        # TODO: Should this use Quantities?
        self.map_to_axis(lambda x: x + offset)

    def map_to_axis(self, fun):
        """
        Maps a function to the given axis. This can be used for non-linear
        corrections of the axis.

        Parameters
        ----------
        fun: Function from float to float
            The function to apply to the wavelengths.
        """
        self.freq_axis = map(fun, self.freq_axis)

    def gaussian_fit(self, guess, *guesses, **kwargs):
        """
        Fits a gaussian distribution to the data, and returns a fit whose
        parameters - amplitude, mean and standard deviation, among others,
        can be called.
        Parameters
        ----------
        guess: tuple of three floats
            The best guess for the first component of the gaussian fit. The
            syntax is (amp_guess, mean_guess, stddev_guess).
        *guesses: additional tuples of three ints
            Additional lines can be fitted by adding more tuples
        **kwargs: dict
            Additional keyword arguments are passed on to the fitter
        """
        g_init = models.Gaussian1D(amplitude=guess[0], mean=guess[1],
                                   stddev=guess[2])
        for (amp, mean, sd) in guesses:
            g_mod = models.Gaussian1D(amplitude=amp, mean=mean, stddev=sd)
            g_init = g_init + g_mod
        fitter = fitting.LevMarLSQFitter()
        return fitter(g_init, self.freq_axis, self.data, **kwargs)
