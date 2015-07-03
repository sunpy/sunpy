# -*- coding: utf-8 -*-
"""
Test for Spectrum
"""
import sunpy.spectra.spectrum as s
from astropy.modeling import models, fitting
import astropy.units as u
import astropy.nddata as ndd
import numpy as np
import pytest


def test_shift_axis():
    axis = np.array([1., 2., 3., 4., 5., 6.])
    spec = s.Spectrum(np.array([0, 0, 0, 0, 0, 0]), axis, u.m)
    spec.shift_axis(0. * u.m)
    assert np.allclose(spec.axis, [1., 2., 3., 4., 5., 6.])
    spec.shift_axis(50 * u.cm)
    assert np.allclose(spec.axis, [1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
    spec.shift_axis(-0.7)
    assert np.allclose(spec.axis, [0.8, 1.8, 2.8, 3.8, 4.8, 5.8])


def test_map_to_axis():
    axis = np.array([1., 2., 3., 4., 5., 6.])
    spec = s.Spectrum(np.array([0, 0, 0, 0, 0, 0]), axis, u.m)
    fun = lambda x: (x - 1 * u.m) ** 2
    spec.map_to_axis(fun)
    assert np.allclose(spec.axis, [0, 1, 4, 9, 16, 25])


def test_gaussian_fit():
    np.random.seed(1)
    axis = np.linspace(-10, 10, 500)
    g = models.Gaussian1D(amplitude=1.5, mean=5, stddev=1.2)
    data = g(axis) + np.random.normal(0, 0.2, 500)
    spec = s.Spectrum(data, axis, u.Angstrom)
    fit = spec.gaussian_fit((1.45, 5.05, 0.3))
    np.allclose(fit.parameters, [1.5, 5, 1.2], 0.05)


def test_gaussian_fit_2():
    np.random.seed(1)
    axis = np.linspace(-10, 10, 500)
    g1 = models.Gaussian1D(amplitude=1.5, mean=5, stddev=1.2)
    g2 = models.Gaussian1D(amplitude=1.0, mean=3, stddev=0.1)
    data = g1(axis) + g2(axis) + np.random.normal(0, 0.2, 500)
    spec = s.Spectrum(data, axis, u.Angstrom)
    fit = spec.gaussian_fit((1.45, 5.05, 0.3), (1, 3, 0.1))
    np.allclose(fit.parameters, [1.5, 5, 1.2, 1.1, 3, 0.1], 0.05)
