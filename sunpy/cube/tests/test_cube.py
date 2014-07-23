# -*- coding: utf-8 -*-
'''
Tests for Cube
'''
from __future__ import absolute_import
import sunpy.cube.cube as c
from sunpy.map.mapbase import GenericMap
from sunpy.spectra.spectrum import Spectrum
from sunpy.spectra.spectrogram import Spectrogram
import numpy as np
from astropy.wcs import WCS
import pytest
import astropy.units as u


# sample data for tests
# TODO: use a fixture reading from a test file. file TBD.
ht = {'CTYPE1': 'HPLT-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5,
      'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2,
      'CTYPE3': 'TIME    ', 'CUNIT3': 'min', 'CDELT3': 0.4}
wt = WCS(header=ht, naxis=3)
data = np.array([[[1,2,3,4], [2,4,5,3], [0,-1,2,3]],
                 [[2,4,5,1], [10,5,2,2], [10,3,3,0]]])
cube = c.Cube(data, wt)

hm = {'CTYPE1': 'HPLT-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5,
      'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2,
      'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4}
wm = WCS(header=hm, naxis=3)
cubem = c.Cube(data, wm)


def test_orient_with_time():
    newdata, newwcs = c._orient(data, wt)
    assert newwcs.wcs.axis_types[-2] == 3000  # code for a spectral dimension
    assert newwcs.wcs.axis_types[-1] == 0  # code for an unknown axis - time
    assert newwcs.naxis == 4
    assert newdata.shape == (2, 3, 4)  # the time dimension should be first
    with pytest.raises(ValueError):
        c._orient(np.zeros((1, 2)), wt)
    with pytest.raises(ValueError):
        c._orient(np.zeros((1, 2, 3, 4)), wt)
    with pytest.raises(ValueError):
        c._orient(data, WCS(naxis=2))


def test_orient_no_time():
    newdata, newwcs = c._orient(data, wm)
    assert newwcs.wcs.axis_types[-1] == 3000
    assert newdata.shape == (3, 2, 4)


def test_slice_to_map_with_time():
    with pytest.raises(c.CubeError):
        cube.slice_to_map(0)
    with pytest.raises(c.CubeError):
        cube.slice_to_map((0, 3))


def test_slice_to_map_no_time():
    m0 = cubem.slice_to_map(0)
    m1 = cubem.slice_to_map((0, 3))
    assert np.all(m0.data == cubem.data[0])
    assert np.all(m1.data == cubem.data.sum(0))


def test_choose_wavelength_slice_with_time():
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


def test_choose_wavelength_no_time():
    ius = cubem._choose_wavelength_slice(-1)  # integer, under range slice
    iis = cubem._choose_wavelength_slice(1)  # integer, in range slice
    ios = cubem._choose_wavelength_slice(11)  # integer, over range slice

    qus = cubem._choose_wavelength_slice(-1 * u.Angstrom)  # quantity, under
    qis = cubem._choose_wavelength_slice(0.4 * u.Angstrom)  # quantity, in
    qos = cubem._choose_wavelength_slice(8 * u.Angstrom)  # quantity, over

    f = cubem._choose_wavelength_slice(0.4)  # no units given

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


def test_select_order():
    lists = [['TIME', 'WAVE', 'HPLT-TAN', 'HPLN-TAN'],
             ['WAVE', 'HPLT-TAN', 'UTC', 'HPLN-TAN'],
             ['HPLT-TAN', 'TIME', 'HPLN-TAN'],
             ['HPLT-TAN', 'DEC--TAN', 'WAVE'],
             [],
             ['UTC', 'TIME', 'WAVE', 'HPLT-TAN']]

    results = [[0, 1, 2, 3],
               [2, 0, 1, 3],
               [1, 0, 2],  # Second order is initial order
               [2, 0, 1],
               [],
               [1, 0, 2, 3]]

    for (l, r) in zip(lists, results):
        assert c._select_order(l) == r


def test_slicing_first_axis():
    # lambda-x-y slices
    s1 = cubem[1]
    s2 = cubem[0:2]
    s3 = cubem[:]

    # time-lambda-y slices
    s4 = cube[1]
    s5 = cube[0:2]
    s6 = cube[:]

    assert isinstance(s1, GenericMap)
    assert isinstance(s2, c.Cube)
    assert isinstance(s3, c.Cube)
    assert isinstance(s4, Spectrum)
    assert isinstance(s5, c.Cube)
    assert isinstance(s6, c.Cube)
    with pytest.raises(IndexError):
        cubem[None]
    with pytest.raises(IndexError):
        cube[None]

    