# Author: Tomas Meszaros <exo@tty.sk>

import astropy.units as u
from sunpy.image.rescale import reshape_image_to_4d_superpixel
import pytest
import os
import numpy as np
import sunpy.data.test

@pytest.fixture
def aia171_test_map():
    testpath = sunpy.data.test.rootdir
    return sunpy.map.Map(os.path.join(testpath, 'aia_171_level1.fits'))

@pytest.fixture
def shape(aia171_test_map):
    return np.array(aia171_test_map.data.shape)

def resample_meta(dimensions, method, center, minusone):
    map_resampled = aia171_test_map().resample(dimensions)
    return tuple(map_resampled.data.shape)

def resample_method(method):
    assert resample_meta((512, 512) * u.pix, method, False, False) == (512, 512)
    assert resample_meta((2056, 2056) * u.pix, method, False, False) == (2056, 2056)
    assert resample_meta((512, 512) * u.pix, method, False, True) == (512, 512)
    assert resample_meta((2056, 2056) * u.pix, method, False, True) == (2056, 2056)
    assert resample_meta((512, 512) * u.pix, method, True, False) == (512, 512)
    assert resample_meta((2056, 2056) * u.pix, method, True, False) == (2056, 2056)
    assert resample_meta((512, 512) * u.pix, method, True, True) == (512, 512)
    assert resample_meta((2056, 2056) * u.pix, method, True, True) == (2056, 2056)

def test_resample_neighbor():
    resample_method('neighbor')

def test_resample_nearest():
    resample_method('nearest')

def test_resample_linear():
    resample_method('linear')

def test_resample_spline():
    resample_method('spline')

def test_reshape(aia171_test_map, shape):
    for factor in np.arange(2, 10):
        if np.all((shape % factor) == 0):
            # it is a factor therefore should work
            reshape_image_to_4d_superpixel(aia171_test_map.data, shape/factor)
        else:
            # it is not a factor so function should raise an error
            with pytest.raises(ValueError):
                reshape_image_to_4d_superpixel(aia171_test_map.data, shape/factor)

