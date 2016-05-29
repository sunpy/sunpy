from __future__ import absolute_import, division, print_function

from sunpy.image import slit
import sunpy.map
import sunpy.data.test
import astropy.units as u
import numpy as np
import os
import pytest
from glob import glob

@pytest.fixture
def eit_test_cube():
    testpath = sunpy.data.test.rootdir
    eit_dir = os.path.join(testpath, 'EIT')
    return sunpy.map.Map(eit_dir, cube=True)

@pytest.fixture
def numbered_array():
    test_numb = np.linspace(1,10,10)*np.ones((10,10))
    return test_numb

@pytest.fixture
def e_list():
    alist=[]
    for i in range(1,11):
        alist.append([i,i])
    return alist


def test_slit():
    a_cube = eit_test_cube():
    test = slit.slit(a_cube, [0, 0]*u.arcsec, [0, 50]*u.arcsec, 0,0,0)
    assert test == a_cube[:, 0, 0:50]

def test_slit_count():
    test_obj = slit.slit_count(0, 0, 928, 0, 0, 1, 1 ) 
    assert isinstance(test_obj, np.ndarray)
    assert slit.slit_count(0, 0, 9, 9, 0, 0, 0) == e_list() 

def test_get_pixels_on_line():
    test_obj = slit.get_pixels_on_line(0, 928, 0, 0) 
    assert isinstance(test_obj, np.ndarray)
    assert slit.get_pixels_on_line(0, 0, 10, 10) == np.linspace(1,10,10)

