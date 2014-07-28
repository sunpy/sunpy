# -*- coding: utf-8 -*-
'''Tests for wcs_util'''
from astropy.wcs import WCS
from sunpy.wcs import wcs_util as wu
import numpy as np


def test_reindex_wcs():
    emptywcs = WCS(naxis=3)
    emptyreindexed = wu.reindex_wcs(emptywcs, np.array([0, 1, 2]))
    assert emptyreindexed.get_axis_types() == emptywcs.get_axis_types()
    h = {'CTYPE1': 'HPLN-TAN', 'CUNIT1': 'deg',
         'CTYPE2': 'HPLT-TAN', 'CUNIT2': 'deg'}
    testwcs = WCS(header=h, naxis=2)
    nochange = wu.reindex_wcs(testwcs, np.array([0, 1]))
    swapped = wu.reindex_wcs(testwcs, np.array([1, 0]))
    assert (testwcs.wcs.ctype[0] == nochange.wcs.ctype[0] and
            testwcs.wcs.ctype[1] == nochange.wcs.ctype[1])
    assert (testwcs.wcs.ctype[0] == swapped.wcs.ctype[1] and
            testwcs.wcs.ctype[1] == swapped.wcs.ctype[0])


def test_add_celestial_axis():
    h = {'CTYPE1': 'TIME    ', 'CUNIT1': 'min', 'CDELT1': 0.5,
         'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2,
         'CTYPE3': 'HPLT-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4}
    w = WCS(header=h, naxis=3, _do_set=False)
    w.wcs.pc = [[1, 2, 3], [3, 0, 1.2], [1, 1, 1]]
    nw = wu.add_celestial_axis(w)
    assert nw.wcs.ctype[3] == 'HPLN-TAN'
    assert nw.wcs.cunit[3] == 'deg'
    assert nw.wcs.cdelt[3] == 1
    assert np.all(nw.wcs.pc == [[1, 2, 3, 0],
                                [3, 0, 1.2, 0],
                                [1, 1, 1, 0],
                                [0, 0, 0, 1]])
