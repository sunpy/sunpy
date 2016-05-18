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

@pytest.fixtre
def numbered_array():
    test_numb = np.linspace(1,10,10)*np.ones((10,10))
    return test_numb

def test_slit_count():
    



def test_slit_count():
    assert slit_count(0, 0, 928, 0, 0, 1, 1 ) isinstance numpy.ndarray
  
    

def test_get_pixels_on_line():
    assert get_pixels_on_line(0, 928, 0, 0) isinstance numpy.ndarray
    assert get_pixels_on_line(0, 0, 10, 10) == np.linspace(1,10,10)




