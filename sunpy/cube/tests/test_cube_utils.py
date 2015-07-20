# -*- coding: utf-8 -*-
import pytest
from sunpy.cube import cube_utils as cu
import numpy as np
from sunpy.wcs.wcs import WCS
from astropy import units as u

ht = {'CTYPE1': 'HPLT-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5, 'CRPIX1': 0, 'CRVAL1': 0,
      'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2, 'CRPIX2': 0, 'CRVAL2': 0,
      'CTYPE3': 'TIME    ', 'CUNIT3': 'min', 'CDELT3': 0.4, 'CRPIX3': 0, 'CRVAL3': 0}
wt = WCS(header=ht, naxis=3)
hm = {'CTYPE1': 'HPLT-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5, 'CRPIX1': 2, 'CRVAL1': 0.5,
      'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2, 'CRPIX2': 0, 'CRVAL2': 10,
      'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1}
wm = WCS(header=hm, naxis=3)
data = np.array([[[1,2,3,4], [2,4,5,3], [0,-1,2,3]],
                 [[2,4,5,1], [10,5,2,2], [10,3,3,0]]])


def test_orient_with_time():
    newdata, newwcs = cu.orient(data, wt)
    assert newwcs.wcs.axis_types[-2] == 3000  # code for a spectral dimension
    assert newwcs.wcs.axis_types[-1] == 0  # code for an unknown axis - time
    assert newwcs.naxis == 4
    assert newdata.shape == (2, 3, 4)  # the time dimension should be first
    with pytest.raises(ValueError):
        cu.orient(np.zeros((1, 2)), wt)
    with pytest.raises(ValueError):
        cu.orient(np.zeros((1, 2, 3, 4)), wt)
    with pytest.raises(ValueError):
        cu.orient(data, WCS(naxis=2))


def test_orient_no_time():
    newdata, newwcs = cu.orient(data, wm)
    assert newwcs.wcs.axis_types[-1] == 3000
    assert newdata.shape == (3, 2, 4)


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
        assert cu.select_order(l) == r


def test_iter_isinstance():
    obj = (1, 'x', 2.5)
    assert cu.iter_isinstance(obj, (int, str, float))
    assert cu.iter_isinstance(obj, ((int, float), (str, int), float))
    assert not cu.iter_isinstance(obj, (int, str))
    assert not cu.iter_isinstance(obj, (int, str, float, int))
    assert not cu.iter_isinstance(1, (int))  # only works for tuples
    assert not cu.iter_isinstance(int, (float))
    assert cu.iter_isinstance(obj, (int, str), (int, str, float), (float, int))


def test_convert_point():
    assert cu.convert_point(10.0, u.Angstrom, wm, 1) == 0
    assert cu.convert_point(10.2, u.Angstrom, wm, 1) == 1
    assert cu.convert_point(9.6, u.Angstrom, wm, 1) == -2
    assert cu.convert_point(10.3, u.Angstrom, wm, 1) == 2
    assert cu.convert_point(0.001, u.mm, wm, 1) == 49950

    assert cu.convert_point(0, u.min, wt, 0) == 0
    assert cu.convert_point(3.1, u.min, wt, 0) == 8
    assert cu.convert_point(-2.4, u.min, wt, 0) == -6
    assert cu.convert_point(0, u.s, wt, 0) == 0
    assert cu.convert_point(24, u.s, wt, 0) == 1
    assert cu.convert_point(-72, u.s, wt, 0) == -3

    assert cu.convert_point(0.2, u.Angstrom, wt, 1) == 1
    assert cu.convert_point(12, None, wt, 0) == 12
    assert cu.convert_point(15.7, None, wm, 3) == 15
    assert cu.convert_point(4, u.pix, wm, 2) == 4
    assert cu.convert_point(7.2, u.pixel, wt, 1) == 7


def test_convert_slice():
    slices = [(slice(1, None, 2.2), wt, 2),
              (slice(9.6, 10.2 * u.Angstrom, None), wm, 1),
              (slice(9.6 * u.Angstrom, 10.2 * u.Angstrom, None), wm, 1),
              (slice(None, None, 0.8 * u.min), wt, 0),
              (slice(3 * u.deg, 14.5, 2), wm, 2)]

    results = [slice(1, None, 2.2),
               slice(-2, 1, None),
               slice(-2, 1, None),
               slice(None, None, 2),
               slice(7, 30, 4)]

    for i in range(len(slices)):
        assert cu._convert_slice(*slices[i]) == results[i]
    with pytest.raises(cu.CubeError):
        cu._convert_slice(slice(2 * u.m, 3, 4 * u.mm), wm, 0)


def test_pixelize_slice():
    sl = [((slice(1, None, 3), 3.2, slice(None, None, 2)), wt),
          ((0.3 * u.deg, 0.0001 * u.mm, 5400 * u.arcsec), wm),
          ((slice(None, 3 * u.deg, 2), slice(5 * u.nm, 6, 1 * u.nm),
            slice(90 * u.arcmin, 300, 60)), wm)]

    results = [(slice(1, None, 3), 3.2, slice(None, None, 2)),
               (0, 4950, 4),
               (slice(None, 7, 5), slice(200, 250, 50), slice(4, 11, 2))]

    for i in range(len(sl)):
        assert cu.pixelize_slice(*sl[i]) == results[i]
