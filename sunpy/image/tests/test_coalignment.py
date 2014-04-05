# Author: Jack Ireland
#
# Testing functions for a mapcube coalignment functionality.  This
# functionality relies on the scikit-image function "match_template".
#

import numpy as np
from numpy.testing import assert_allclose
from scipy.ndimage.interpolation import shift
from sunpy import AIA_171_IMAGE
from sunpy import map
from sunpy.image.coalignment import parabolic_turning_point, \
repair_nonfinite, default_fmap_function, _lower_clip, _upper_clip, \
calculate_clipping, get_correlation_shifts, find_best_match_location, \
match_template_to_layer, clip_edges, calculate_shift


# Map and template we will use in testing
testmap = map.Map(AIA_171_IMAGE)
test_layer = testmap.data
ny = test_layer.shape[0]
nx = test_layer.shape[1]
test_template = test_layer[1 + ny / 4 : 1 + 3 * ny / 4,
                            2 + nx / 4 : 2 + 3 * nx / 4]

# Used in testing the clipping
clip_test_array = np.asarray([0.2, -0.3, -1.0001])


def test_parabolic_turning_point():
    assert(parabolic_turning_point(np.asarray([6.0, 2.0, 0.0])) == 1.5)


def test_repair_nonfinite():
    a = np.ones((9))
    for i in range(0, 9):
        for non_number in [np.nan, np.inf]:
            a[i] = non_number
            b = a.reshape(3, 3)
            c = repair_nonfinite(b)
            assert(np.isfinite(c).all())


def test_match_template_to_layer():
    result = match_template_to_layer(test_layer, test_template)
    assert(result.shape[0] == 513)
    assert(result.shape[1] == 513)
    assert_allclose(np.max(result), 1.00, rtol=1e-2, atol=0)


def test_get_correlation_shifts():
    test_array = np.zeros((3,3))
    test_array[1, 1] = 1
    test_array[2, 1] = 0.6
    test_array[1, 2] = 0.2
    y_test, x_test = get_correlation_shifts(test_array)
    assert_allclose(y_test, 0.214285714286, rtol=1e-2, atol=0)
    assert_allclose(x_test, 0.0555555555556, rtol=1e-2, atol=0)


def test_find_best_match_location():
    result = match_template_to_layer(test_layer, test_template)
    y_test, x_test = find_best_match_location(result)
    assert_allclose(y_test, 257.0, rtol=1e-3, atol=0)
    assert_allclose(x_test, 258.0, rtol=1e-3, atol=0)

def test_lower_clip():
    assert(_lower_clip(clip_test_array) == 2.0)


def test_upper_clip():
    assert(_upper_clip(clip_test_array) == 1.0)


def test_calculate_clipping():
    answer = calculate_clipping(clip_test_array, clip_test_array)
    assert(answer == ([2.0, 1.0], [2.0, 1.0]))


def test_clip_edges():
    a = np.zeros(shape=(341, 156))
    yclip = [4, 0]
    xclip = [1, 2]
    new_a = clip_edges(a, yclip, xclip)
    assert(a.shape[0] - (yclip[0] + yclip[1]) == 337)
    assert(a.shape[1] - (xclip[0] + xclip[1]) == 153)


def test_calculate_shift():
    result = calculate_shift(test_layer, test_template)
    assert_allclose(result[0], 257.0,  rtol=1e-3, atol=0)
    assert_allclose(result[1], 258.0,  rtol=1e-3, atol=0)


def test_default_fmap_function():
    assert(default_fmap_function([1,2,3]).dtype == np.float64(1).dtype)


def test_coalign():
    # take the AIA image and shift it
    # Pixel displacements have the y-displacement as the first entry
    pixel_displacements = np.asarray([1.6, 10.1])
    known_displacements = {'x':np.asarray([0.0, pixel_displacements[1] * testmap.scale['x']]), 'y':np.asarray([0.0, pixel_displacements[0] * testmap.scale['y']])}
    # Apply the shift
    d1 = shift(testmap.data, pixel_displacements)
    m1 = map.Map((d1, testmap.meta))
    # Create the mapcube
    mc = map.Map([testmap, m1], cube=True)
    # Do the coalignment
    displacements = mc.coalign(displacements_only=True)
    # Assert
    assert_allclose(displacements['x'], known_displacements['x'], rtol=5e-2, atol=0)
    assert_allclose(displacements['y'], known_displacements['y'], rtol=5e-2, atol=0 )

