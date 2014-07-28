# -*- coding: utf-8 -*-
'''
Tests for Cube
'''
from __future__ import absolute_import
import sunpy.cube.cube as c
import sunpy.cube.cube_utils as cu
from sunpy.map.mapbase import GenericMap
from sunpy.spectra.spectrum import Spectrum
from sunpy.spectra.spectrogram import Spectrogram
from sunpy.lightcurve.lightcurve import LightCurve
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


def test_slice_to_map_with_time():
    with pytest.raises(cu.CubeError):
        cube.slice_to_map(0)
    with pytest.raises(cu.CubeError):
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


def test_freq_axis():
    f1 = cube.freq_axis()
    f2 = cubem.freq_axis()
    # the e-11 are the conversions from angstrom to meters
    assert np.allclose(f1, [0, 2.0e-11, 4.0e-11])
    assert np.allclose(f2, [0, 2.0e-11, 4.0e-11])


def test_time_axis():
    t1 = cube.time_axis()
    assert np.allclose(t1, [0, 0.4])
    with pytest.raises(cu.CubeError):
        cubem.time_axis()


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


def test_slicing_second_axis():
    # lambda-x-y
    slices = [cubem[:, 1], cubem[:, 0:2], cubem[:, :],
              cubem[1, 1], cubem[1, 0:2], cubem[1, :],
              # time-lambda-y
              cube[:, 1], cube[:, 0:2], cube[:, :],
              cube[1, 1], cube[1, 0:2], cube[1, :]]

    types = [np.ndarray, c.Cube, c.Cube, np.ndarray, GenericMap, GenericMap,
             LightCurve, c.Cube, c.Cube, np.ndarray, Spectrum, Spectrum]
    for (s, t) in zip(slices, types):
        assert isinstance(s, t)


def test_slicing_third_axis():
    slices = [cubem[:, :, 1], cubem[:, :, 0:2], cubem[:, :, :],
              cubem[:, 1, 1], cubem[:, 1, 0:2], cubem[:, 1, :],
              cubem[1, :, 1], cubem[1, :, 0:2], cubem[1, :, :],
              cubem[1, 1, 1], cubem[1, 1, 0:2], cubem[1, 1, :],
              # time-lambda-y
              cube[:, :, 1], cube[:, :, 0:2], cube[:, :, :],
              cube[:, 1, 1], cube[:, 1, 0:2], cube[:, 1, :],
              cube[1, :, 1], cube[1, :, 0:2], cube[1, :, :],
              cube[1, 1, 1], cube[1, 1, 0:2], cube[1, 1, :]]

    types = [np.ndarray, c.Cube, c.Cube,
             Spectrum, np.ndarray, np.ndarray,
             np.ndarray, GenericMap, GenericMap,
             int, np.ndarray, np.ndarray,

             Spectrogram, c.Cube, c.Cube,
             LightCurve, LightCurve, LightCurve,
             Spectrum, Spectrum, Spectrum,
             int, np.ndarray, np.ndarray]
    for (s, t) in zip(slices, types):
        assert isinstance(s, t)


def test_reduce_dim():
    slices = [slice(s, e, t) for s, e, t in [(None, None, None), (0, 2, None),
                                             (None, None, 2)]]
    assert np.all(cube.data == cube._reduce_dim(0, slices[0]))
    assert np.all(cube.data == cube._reduce_dim(0, slices[1]))
    assert cubem._reduce_dim(0, slices[1]).data.shape == (2, 2, 4)
    assert cubem._reduce_dim(2, slices[2]).data.shape == (3, 2, 2)
    assert cube._reduce_dim(2, slices[2]).axes_wcs.wcs.cdelt[1] == 1
    