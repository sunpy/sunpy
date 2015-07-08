# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>,
#         Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101, W0141
"""
Module containing the definition of the Spectrum object. Refer to the Spectrum
documentation for more details.
"""
from __future__ import absolute_import

import astropy.nddata as ndd
from astropy.modeling import models, fitting
import astropy.units as u
import numpy as np

from matplotlib import pyplot as plt

__all__ = ['Spectrum']


class Spectrum(ndd.NDDataArray):
    """
    Class representing a spectrum.

    Attributes
    ----------
    axis: np.ndarray
        One-dimensional array with the frequency or wavelength values at every
        data point.

    axis_unit: astropy unit
        The unit of the spectral axis. This must be in units of frequency or
        distance.

    data: np.ndarray
        one-dimensional array which the intensity at a particular frequency at
        every data-point.
    """

    def __init__(self, data, axis, axis_unit, **kwargs):
        ndd.NDDataArray.__init__(self, data=data, **kwargs)
        self.axis = axis
        self.axis_unit = axis_unit

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
            lines = axes.plot(self.axis, self, **params)
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
        offset: float or astropy Quantity
            The amount to offset by. If no unit is given the current axis unit
            is used
        """
        if isinstance(offset, u.Quantity):
            self.map_to_axis(lambda x: x + offset)
        else:
            self.map_to_axis(lambda x: x + (offset * self.axis_unit))

    def map_to_axis(self, fun):
        """
        Maps a function to the given axis. This can be used for non-linear
        corrections of the axis.

        Parameters
        ----------
        fun: Function from Quantity to Quantity
            The function to apply to the wavelengths.
        """
        qtys = [tick * self.axis_unit for tick in self.axis]
        newqtys = map(fun, qtys)
        newaxis = [tick.value for tick in newqtys]
        self.axis = newaxis

    def gaussian_fit(self, line_guess, *extra_lines, **kwargs):
        """
        Fits a gaussian distribution to the data, and returns a fit whose
        parameters - amplitude, mean and standard deviation, among others,
        can be called.

        Parameters
        ----------
        line_guess: tuple of three floats
            The best guess for the first component of the gaussian fit. The
            syntax is (amp_guess, mean_guess, stddev_guess).
        *extra_lines: additional tuples of three ints
            Additional lines can be fitted by adding more tuples
        **kwargs: dict
            Additional keyword arguments are passed on to the fitter
        """
        g_init = models.Gaussian1D(amplitude=line_guess[0], mean=line_guess[1],
                                   stddev=line_guess[2])
        for (amp, mean, stddev) in extra_lines:
            g_mod = models.Gaussian1D(amplitude=amp, mean=mean, stddev=stddev)
            g_init = g_init + g_mod
        fitter = fitting.LevMarLSQFitter()
        fit_axis = np.empty()
        fit_data = np.empty()
        x_range = kwargs.get('x_range')
        if x_range is not None:
            arrmin = self.axis.index(x_range[0] / self.axis_unit)
            arrmax = self.axis.index(x_range[1] / self.axis_unit)
            fit_axis = self.axis[arrmin:arrmax]
            fit_data = self.data[arrmin:arrmax]
        else:
            fit_axis = self.axis
            fit_data = self.data
        return fitter(g_init, fit_axis, fit_data, **kwargs)

    def __getitem__(self, item):
        # TODO: use quantities here
        newdata = self.data.__getitem__(item)
        newaxis = self.axis.__getitem__(item)
        return Spectrum(newdata, newaxis, self.axis_unit)
