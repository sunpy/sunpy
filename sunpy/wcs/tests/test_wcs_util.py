# -*- coding: utf-8 -*-
'''Tests for wcs_util'''
from astropy.wcs import WCS
from sunpy.util import wcs_util as wu
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