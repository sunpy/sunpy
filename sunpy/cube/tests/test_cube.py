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
from sunpy.wcs.wcs import WCS
import pytest
import astropy.units as u


# sample data for tests
# TODO: use a fixture reading from a test file. file TBD.
ht = {'CTYPE1': 'HPLT-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5, 'CRPIX1': 0, 'CRVAL1': 0,
      'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2, 'CRPIX2': 0, 'CRVAL2': 0,
      'CTYPE3': 'TIME    ', 'CUNIT3': 'min', 'CDELT3': 0.4, 'CRPIX3': 0, 'CRVAL3': 0}
wt = WCS(header=ht, naxis=3)
data = np.array([[[1,2,3,4], [2,4,5,3], [0,-1,2,3]],
                 [[2,4,5,1], [10,5,2,2], [10,3,3,0]]])
cube = c.Cube(data, wt)

hm = {'CTYPE1': 'HPLT-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5, 'CRPIX1': 2, 'CRVAL1': 0.5,
      'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2, 'CRPIX2': 0, 'CRVAL2': 10,
      'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1}
wm = WCS(header=hm, naxis=3)
cubem = c.Cube(data, wm)

h4 = {'CTYPE1': 'HPLT-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5, 'CRPIX1': 0, 'CRVAL1': 0,
      'CTYPE2': 'WAVE    ', 'CUNIT2': 'nm', 'CDELT2': 0.2, 'CRPIX2': 0, 'CRVAL2': 10,
      'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1,
      'CTYPE4': 'TIME    ', 'CUNIT4': 'min', 'CDELT4': 0.6, 'CRPIX4': 0, 'CRVAL4': 0}
w4 = WCS(header=h4, naxis=4)
data4 = np.array([[[[1,2,3,4], [2,4,5,3], [0,-1,2,3]],
                   [[2,4,5,1], [10,5,2,2], [10,3,3,0]]],

                  [[[4,6,3,5], [1,0,5,3], [-4,8,7,6]],
                   [[7,3,2,6], [1,7,8,7], [2,4,0,1]]]])
hcube = c.Cube(data4, w4)

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
    f1, u1 = cube.freq_axis()
    f2, u2 = cubem.freq_axis()
    # the e-11 are the conversions from angstrom to meters
    assert np.allclose(f1, [0, 2.0e-11, 4.0e-11])
    assert np.allclose(f2, [0, 2.0e-11, 4.0e-11])
    assert u1 == u.m
    assert u2 == u.m


def test_time_axis():
    t1, u1 = cube.time_axis()
    assert np.allclose(t1, [0, 0.4])
    assert u1 == u.min
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


def test_slicing_first_axis_world_coord():
    # lambda-x-y slices
    s1 = cubem[10.2 * u.Angstrom]
    s2 = cubem[10 * u.Angstrom:10.4]
    s3 = cubem[::0.4 * u.Angstrom]

    # time-lambda-y slices
    s4 = cube[0.5 * u.min]
    s5 = cube[0:0.8 * u.min]
    s6 = cube[::0.8 * u.min]

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


def test_slicing_second_axis_world_coord():
    # lambda-x-y
    slices = [cubem[:, -0.5 * u.deg], cubem[:, -0.5:0.5 * u.deg],
              cubem[:, ::1 * u.deg], cubem[10.2 * u.Angstrom, 0 * u.deg],
              cubem[1, -0.5:0.5 * u.deg], cubem[1, ::1 * u.deg],
              # time-lambda-y
              cube[:, 0.2 * u.Angstrom], cube[:, 0:0.4 * u.Angstrom],
              cube[:, 0:0.4 * u.Angstrom:0.4], cube[1, 0.2 * u.Angstrom],
              cube[0.5 * u.min, 0:2], cube[0.5 * u.min, 0:0.4 * u.Angstrom:0.4]
              ]

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


def test_slicing_third_axis_world_coord():
    slices = [cubem[:, :, 0.6 * u.deg], cubem[:, :, 0.2:1 * u.deg],
              cubem[:, :, 0.2:1 * u.deg:0.8],
              cubem[:, -0.5 * u.deg, 0.6 * u.deg], cubem[:, -0.5 * u.deg, 0:2],
              cubem[:, -0.5 * u.deg, 0.2:1 * u.deg:0.8],
              cubem[1.02 * u.nm, :, 0.6 * u.deg], cubem[1, :, 0.2:1 * u.deg],
              cubem[1.02 * u.nm, :, 0.2:1 * u.deg:0.8],
              cubem[1.02 * u.nm, -0.5 * u.deg, 0.6 * u.deg],
              cubem[1, 1, 0.2:1 * u.deg],
              cubem[1.02 * u.nm, -0.5 * u.deg, 0.2:1 * u.deg:0.8],
              # time-lambda-y
              cube[:, :, 0.5 * u.deg], cube[:, :, 0:1 * u.deg],
              cube[:, :, ::1 * u.deg], cube[:, 0.2 * u.Angstrom, 0.5 * u.deg],
              cube[:, 1, 0:1 * u.deg], cube[:, 0.2 * u.Angstrom, ::1 * u.deg],
              cube[0.5 * u.min, :, 0.5 * u.deg], cube[1, :, 0:1 * u.deg],
              cube[0.5 * u.min, :, ::1 * u.deg],
              cube[0.5 * u.min, 0.2 * u.Angstrom, 0.5 * u.deg],
              cube[0.5 * u.min, 1, 0:1 * u.deg], cube[0.5 * u.min, 1, :]]

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


def test_4d_getitem_to_array():
    slices = [hcube[1, 1, 1, 1], hcube[0, 1, 1],
              hcube[2, 0, 1, :], hcube[3, 1, :, 2]]
    assert isinstance(slices[0], int)
    for s in slices[1:]:
        assert isinstance(s, np.ndarray)


def test_4d_getitem_to_array_world_coord():
    slices = [hcube[1, 10.2 * u.nm, 1, 0.6 * u.deg], hcube[3, 1, :, 1 * u.deg],
              hcube[0.7 * u.min, 1, 0.5 * u.deg],
              hcube[1 * u.min, 10 * u.nm, 1, :]]
    assert isinstance(slices[0], int)
    for s in slices[1:]:
        assert isinstance(s, np.ndarray)


def test_4d_getitem_to_map():
    slices = [hcube[2, 0], hcube[1, 1, :], hcube[1, 1, :, 0:2]]
    for s in slices:
        assert isinstance(s, GenericMap)


def test_4d_getitem_to_map_world_coord():
    slices = [hcube[1 * u.min, 0], hcube[0.4 * u.min, 10.2 * u.nm, :],
              hcube[0.4 * u.min, 10.2 * u.nm, :, 0:2]]
    for s in slices:
        assert isinstance(s, GenericMap)


def test_4d_getitem_to_spectrum():
    slices = [hcube[1, :, 1, 2], hcube[3, :, 0], hcube[2, :, 1, 0:2],
              hcube[0, :, :, 0]]
    for s in slices:
        assert isinstance(s, Spectrum)


def test_4d_getitem_to_spectrum_world_coord():
    slices = [hcube[0.4 * u.min, :, 0.7 * u.deg, 2],
              hcube[1.3 * u.min, :, 0.7 * u.deg], hcube[0.8 * u.min, :, 1, :2],
              hcube[0 * u.min, :, :, 0]]
    for s in slices:
        assert isinstance(s, Spectrum)


def test_4d_getitem_to_cube():
    slices = [hcube[2], hcube[1, 0:1], hcube[3, :, 0:2], hcube[0, :, :, 0:2],
              hcube[1:3, 1], hcube[0:2, 0, :], hcube[:, 0, :, 0:2],
              hcube[1:3, :, 1, 1:2], hcube[:, :, 0], hcube[:, :, :, 2]]
    for s in slices:
        assert isinstance(s, c.Cube) and s.data.ndim == 3


def test_4d_getitem_to_cube_world_coord():
    slices = [hcube[1.3 * u.min], hcube[1.3 * u.min, 0:1],
              hcube[1.3 * u.min, :, 0:2], hcube[1.3 * u.min, :, :, 0:2],
              hcube[1:3, 10.2 * u.nm], hcube[0:2, 10.2 * u.nm, :],
              hcube[:, 10.2 * u.nm, :, 0:2],
              hcube[0.4 * u.min:1.2 * u.min:0.8, :, 1, 1:2],
              hcube[:, :, 0.7 * u.deg], hcube[:, :, :, 0.6 * u.deg]]
    for s in slices:
        assert isinstance(s, c.Cube) and s.data.ndim == 3


def test_4d_getitem_to_hypercube():
    slices = [hcube[1:3, :, :, 0:1], hcube[1:, :1, :], hcube[2:, :], hcube[1:]]
    for s in slices:
        assert isinstance(s, c.Cube) and s.data.ndim == 4


def test_4d_getitem_to_hypercube_world_coord():
    slices = [hcube[0.4 * u.min:1.2 * u.min:0.8, :, :, 0:1],
              hcube[1:, 10:10.4 * u.nm:2, :], hcube[2:, 10:10.4 * u.nm:2],
              hcube[0.4 * u.min:1.2 * u.min:0.8]]
    for s in slices:
        assert isinstance(s, c.Cube) and s.data.ndim == 4


def test_4d_getitem_to_spectrogram():
    s = hcube[2:, :, 1, 2]
    assert isinstance(s, Spectrogram)


def test_4d_getitem_to_spectrogram_world_coord():
    s = hcube[0.4 * u.min:1.2:0.8, 10:10.4 * u.nm:2, 1, 0.6 * u.deg]
    assert isinstance(s, Spectrogram)


def test_4d_getitem_to_lightcurve():
    slices = [hcube[:, 0, 0, 0], hcube[:, 1, 1, :], hcube[:, 1, 0],
              hcube[:, 1, :, 2]]
    for s in slices:
        assert isinstance(s, LightCurve)


def test_4d_getitem_to_lightcurve_world_coord():
    slices = [hcube[:, 0, 0, 0.6 * u.deg], hcube[:, 1, :, 1 * u.deg],
              hcube[0.4 * u.min:1.2:0.8, 1, 0.5 * u.deg, :],
              hcube[0.4 * u.min:1.2:0.8, 1, 0.6 * u.deg]]
    for s in slices:
        assert isinstance(s, LightCurve)


def test_reduce_dim():
    slices = [slice(s, e, t) for s, e, t in [(None, None, None), (0, 2, None),
                                             (None, None, 2)]]
    assert np.all(cube.data == cu.reduce_dim(cube, 0, slices[0]))
    assert np.all(cube.data == cu.reduce_dim(cube, 0, slices[1]))
    assert cu.reduce_dim(cubem, 0, slices[1]).data.shape == (2, 2, 4)
    assert cu.reduce_dim(cubem, 2, slices[2]).data.shape == (3, 2, 2)
    assert cu.reduce_dim(cube, 2, slices[2]).axes_wcs.wcs.cdelt[1] == 1
