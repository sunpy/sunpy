# -*- coding: utf-8 -*-
'''
Tests for SpectralCube
'''
from __future__ import absolute_import
import sunpy.cube.spectral_cube as sc
import numpy as np
from astropy.wcs import WCS
import pytest
import astropy.units as u


# sample data for tests
# TODO: use a fixture reading from a test file. file TBD.
h = {'CTYPE1': 'HPLN-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5,
     'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2,
     'CTYPE3': 'HPLT-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4}
w = WCS(header=h, naxis=3)
data = np.array([[[1,2,3,4], [2,4,5,3], [0,-1,2,3]],
                 [[2,4,5,1], [10,5,2,2], [10,3,3,0]]])
cube = sc.SpectralCube(data, w)


def test_orient():
    newdata, newwcs = sc._orient(data, w)
    assert newwcs.wcs.axis_types[2] == 3000  # code for a spectral dimension
    assert newdata.shape == (3, 2, 4)  # the spectral dimension should be first
    with pytest.raises(ValueError):
        sc._orient(np.zeros((1, 2)), w)
        sc._orient(np.zeros((1, 2, 3, 4)), w)
        sc._orient(data, WCS(naxis=2))


def test_slice_to_map():
    m0 = cube.slice_to_map(0)
    m1 = cube.slice_to_map((0, 3))
    assert np.all(m0.data == cube.data[0])
    assert np.all(m1.data == cube.data.sum(0))


def test_choose_wavelength_slice():
    ius = cube._choose_wavelength_slice(-1)  # integer, under range slice
    iis = cube._choose_wavelength_slice(1)  # integer, in range slice
    ios = cube._choose_wavelength_slice(11)  # integer, over range slice

    qus = cube._choose_wavelength_slice(-1 * u.Angstrom)  # quantity, under
    qis = cube._choose_wavelength_slice(0.5 * u.Angstrom)  # quantity, in
    qos = cube._choose_wavelength_slice(8 * u.Angstrom)  # quantity, over

    f = cube._choose_wavelength_slice(0.4)  # no units given

    assert ius is None
    assert np.all(iis == [[2, 4, 5, 3], [10, 5, 2, 2]])
    assert ios is None

    assert qus is None
    assert np.all(qis == [[0, -1, 2, 3], [10, 3, 3, 0]])
    assert qos is None

    assert f is None


def test_choose_x_slice():
    ius = cube._choose_x_slice(-1)  # integer, under range slice
    iis = cube._choose_x_slice(1)  # integer, in range slice
    ios = cube._choose_x_slice(11)  # integer, over range slice

    qus = cube._choose_x_slice(-1 * u.deg)  # quantity, under
    qis = cube._choose_x_slice(1.5 * u.deg)  # quantity, in
    qos = cube._choose_x_slice(8 * u.deg)  # quantity, over

    f = cube._choose_x_slice(0.4)  # no units given

    assert ius is None
    assert np.all(iis == [[2, 4, -1], [4, 5, 3]])
    assert ios is None

    assert qus is None
    assert np.all(qis == [[4, 3, 3], [1, 2, 0]])
    assert qos is None

    assert f is None
