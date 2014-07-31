# -*- coding: utf-8 -*-
import pytest
from sunpy.cube import cube_utils as cu
import numpy as np
from sunpy.wcs.wcs import WCS

ht = {'CTYPE1': 'HPLT-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5,
      'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2,
      'CTYPE3': 'TIME    ', 'CUNIT3': 'min', 'CDELT3': 0.4}
wt = WCS(header=ht, naxis=3)
hm = {'CTYPE1': 'HPLT-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5,
      'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2,
      'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4}
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
    assert cu.iter_isinstance(obj, int, str, float)
    assert cu.iter_isinstance(obj, (int, float), (str, int), float)
    assert not cu.iter_isinstance(obj, int, str)
    assert not cu.iter_isinstance(obj, int, str, float, int)
    assert not cu.iter_isinstance(1, int)  # only works for tuples
    assert not cu.iter_isinstance(int, float)
