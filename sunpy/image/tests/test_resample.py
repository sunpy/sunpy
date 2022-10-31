
import numpy as np
import pytest

import astropy.units as u

from sunpy.image.resample import reshape_image_to_4d_superpixel


@pytest.fixture
def shape(aia171_test_map):
    return np.array(aia171_test_map.data.shape)


def resample_meta(aia171_test_map, dimensions, method, center, minusone):
    map_resampled = aia171_test_map.resample(dimensions, method=method)
    return tuple(map_resampled.data.shape)


def resample_method(aia171_test_map, method):
    for shape in [(64, 64), (256, 256)]:
        for center in [False, True]:
            for minusone in [False, True]:
                assert resample_meta(aia171_test_map, shape * u.pix, method, center, minusone) == shape


def test_resample_nearest(aia171_test_map):
    resample_method(aia171_test_map, 'nearest')


def test_resample_linear(aia171_test_map):
    resample_method(aia171_test_map, 'linear')


def test_resample_spline(aia171_test_map):
    resample_method(aia171_test_map, 'spline')


def test_reshape(aia171_test_map, shape):

    def _n(a, b, c):
        return int(np.floor((a-b)/c))

    # Dimension divides the array shape exactly with no remainder
    im = reshape_image_to_4d_superpixel(aia171_test_map.data, (2, 2), (0, 0))
    assert im.shape == (shape[0]/2, 2, shape[1]/2, 2)
    # Dimension divides the array shape exactly with remainder
    im = reshape_image_to_4d_superpixel(aia171_test_map.data, (7, 5), (0, 0))
    assert im.shape == (int(shape[0]/7), 7, int(shape[1]/5), 5)
    # Dimension divides the array shape exactly with no remainder, and there is
    # an offset
    im = reshape_image_to_4d_superpixel(aia171_test_map.data, (2, 2), (1, 1))
    assert im.shape == (int(shape[0]/2) - 1, 2, int(shape[1]/2) - 1, 2)
    # Dimension divides the array shape exactly with remainder, and there is
    # an offset
    d = (9, 7)
    o = (1, 4)
    im = reshape_image_to_4d_superpixel(aia171_test_map.data, d, o)
    assert im.shape == (_n(shape[0], o[0], d[0]), d[0],
                        _n(shape[1], o[1], d[1]), d[1])
    im = reshape_image_to_4d_superpixel(aia171_test_map.data, d, o)
    assert im.shape == (_n(shape[0], o[0], d[0]), d[0],
                        _n(shape[1], o[1], d[1]), d[1])

    d = (9, 7)
    o = (5, 4)
    im = reshape_image_to_4d_superpixel(aia171_test_map.data, d, o)
    assert im.shape == (_n(shape[0], o[0], d[0]), d[0],
                        _n(shape[1], o[1], d[1]), d[1])

    d = (9, 7)
    o = (4, 4)
    im = reshape_image_to_4d_superpixel(aia171_test_map.data, d, o)
    assert im.shape == (_n(shape[0], o[0], d[0]), d[0],
                        _n(shape[1], o[1], d[1]), d[1])
