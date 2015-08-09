# -*- coding: utf-8 -*-
"""
Test for Spectrum
"""
import sunpy.spectra.spectrum as s
import numpy as np
import pytest

axis = np.array([1., 2., 3., 4., 5., 6.])
spec = s.Spectrum(np.array([0, 0, 0, 0, 0, 0]), axis)


def test_shift_axis():
    spec.shift_axis(0.)
    assert np.allclose(spec.freq_axis, [1., 2., 3., 4., 5., 6.])
