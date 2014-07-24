# Author: Jack Ireland
#
# Testing functions for a mapcube coalignment functionality.  This
# functionality relies on the scikit-image function "match_template".
#

import numpy as np
from astropy import units as u
from numpy.testing import assert_allclose, assert_array_almost_equal
from scipy.ndimage.interpolation import shift
from sunpy import AIA_171_IMAGE
from sunpy import map
from sunpy.image.coalignment import parabolic_turning_point, \
repair_image_nonfinite, _default_fmap_function, _lower_clip, _upper_clip, \
calculate_clipping, get_correlation_shifts, find_best_match_location, \
match_template_to_layer, clip_edges, calculate_shift, \
mapcube_coalign_by_match_template

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


def test_repair_image_nonfinite():
    for i in range(0, 9):
        for non_number in [np.nan, np.inf]:
            a = np.ones((9))
            a[i] = non_number
            b = a.reshape(3, 3)
            c = repair_image_nonfinite(b)
            assert(np.isfinite(c).all())


def test_match_template_to_layer():
    result = match_template_to_layer(test_layer, test_template)
    assert(result.shape[0] == 513)
    assert(result.shape[1] == 513)
    assert_allclose(np.max(result), 1.00, rtol=1e-2, atol=0)


def test_get_correlation_shifts():
    # Input array is 3 by 3, the most common case
    test_array = np.zeros((3, 3))
    test_array[1, 1] = 1
    test_array[2, 1] = 0.6
    test_array[1, 2] = 0.2
    y_test, x_test = get_correlation_shifts(test_array)
    assert_allclose(y_test, 0.214285714286, rtol=1e-2, atol=0)
    assert_allclose(x_test, 0.0555555555556, rtol=1e-2, atol=0)

    # Input array is smaller in one direction than the other.
    test_array = np.zeros((2, 2))
    test_array[0, 0] = 0.1
    test_array[0, 1] = 0.2
    test_array[1, 0] = 0.4
    test_array[1, 1] = 0.3
    y_test, x_test = get_correlation_shifts(test_array)
    assert_allclose(y_test, 1.0, rtol=1e-2, atol=0)
    assert_allclose(x_test, 0.0, rtol=1e-2, atol=0)

    # Input array is too big in either direction
    test_array = np.zeros((4, 3))
    y_test, x_test = get_correlation_shifts(test_array)
    assert(y_test == None)
    assert(x_test == None)
    test_array = np.zeros((3, 4))
    y_test, x_test = get_correlation_shifts(test_array)
    assert(y_test == None)
    assert(x_test == None)

def test_find_best_match_location():
    result = match_template_to_layer(test_layer, test_template)
    y_test, x_test = find_best_match_location(result)
    assert_allclose(y_test, 257.0, rtol=1e-3, atol=0)
    assert_allclose(x_test, 258.0, rtol=1e-3, atol=0)

def test_lower_clip():
    assert(_lower_clip(clip_test_array) == 2.0)
    # No element is less than zero
    test_array = np.asarray([1.1, 0.1, 3.0])
    assert(_lower_clip(test_array) == 0)


def test_upper_clip():
    assert(_upper_clip(clip_test_array) == 1.0)
    # No element is greater than zero
    test_array = np.asarray([-1.1, -0.1, -3.0])
    assert(_upper_clip(test_array) == 0)


def test_calculate_clipping():
    answer = calculate_clipping(clip_test_array *u.pix, clip_test_array *u.pix)
    assert_array_almost_equal(answer, ([2.0, 1.0]*u.pix, [2.0, 1.0]*u.pix))


def test_clip_edges():
    a = np.zeros(shape=(341, 156))
    yclip = [4, 0] * u.pix
    xclip = [1, 2] * u.pix
    new_a = clip_edges(a, yclip, xclip)
    assert(a.shape[0] - (yclip[0].value + yclip[1].value) == 337)
    assert(a.shape[1] - (xclip[0].value + xclip[1].value) == 153)


def test_calculate_shift():
    result = calculate_shift(test_layer, test_template)
    assert_allclose(result[0], 257.0,  rtol=1e-3, atol=0)
    assert_allclose(result[1], 258.0,  rtol=1e-3, atol=0)


def test__default_fmap_function():
    assert(_default_fmap_function([1,2,3]).dtype == np.float64(1).dtype)


def test_mapcube_coalign_by_match_template():
    # take the AIA image and shift it
    # Pixel displacements have the y-displacement as the first entry
    pixel_displacements = np.asarray([1.6, 10.1])
    known_displacements = {'x':np.asarray([0.0, pixel_displacements[1] * testmap.scale['x']]), 'y':np.asarray([0.0, pixel_displacements[0] * testmap.scale['y']])}

    # Create a map that has been shifted a known amount.
    d1 = shift(testmap.data, pixel_displacements)
    m1 = map.Map((d1, testmap.meta))

    # Create the mapcube
    mc = map.Map([testmap, m1], cube=True)

    # Test to see if the code can recover the displacements. Do the coalignment
    # using the "return_displacements_only" option
    test_displacements = mapcube_coalign_by_match_template(mc, return_displacements_only=True)
    # Assert
    assert_allclose(test_displacements['x'], known_displacements['x'], rtol=5e-2, atol=0)
    assert_allclose(test_displacements['y'], known_displacements['y'], rtol=5e-2, atol=0 )

    # Test setting the template as a ndarray
    template_ndarray = testmap.data[ny / 4: 3 * ny / 4, nx / 4: 3 * nx / 4]
    test_displacements = mapcube_coalign_by_match_template(mc, template=template_ndarray, return_displacements_only=True)
    # Assert
    assert_allclose(test_displacements['x'], known_displacements['x'], rtol=5e-2, atol=0)
    assert_allclose(test_displacements['y'], known_displacements['y'], rtol=5e-2, atol=0 )

    # Test setting the template as GenericMap
    submap = testmap.submap([nx / 4, 3 * nx / 4], [ny / 4, 3 * ny / 4], units='pixels')
    test_displacements = mapcube_coalign_by_match_template(mc, template=submap, return_displacements_only=True)
    # Assert
    assert_allclose(test_displacements['x'], known_displacements['x'], rtol=5e-2, atol=0)
    assert_allclose(test_displacements['y'], known_displacements['y'], rtol=5e-2, atol=0 )

    # Test setting the template as something other than a ndarray and a
    # GenericMap.  This should throw a ValueError.
    try:
        test_displacements = mapcube_coalign_by_match_template(mc, template='broken')
    except ValueError:
        pass

    # Test passing in displacements
    test_apply_displacements = {'x':-test_displacements['x'], 'y':-test_displacements['y']}
    test_displacements = mapcube_coalign_by_match_template(mc,
                                                           apply_displacements=test_apply_displacements,
                                                           return_displacements_only=True)
    assert_allclose(test_displacements['x'], test_apply_displacements['x'], rtol=5e-2, atol=0)
    assert_allclose(test_displacements['y'], test_apply_displacements['y'], rtol=5e-2, atol=0)

    # Test returning using the "with_displacements" option
    test_output = mapcube_coalign_by_match_template(mc, with_displacements=True)
    # Assert
    assert(isinstance(test_output[0], map.MapCube))
    assert_allclose(test_output[1]['x'], known_displacements['x'], rtol=5e-2, atol=0)
    assert_allclose(test_output[1]['y'], known_displacements['y'], rtol=5e-2, atol=0 )

    # Test returning with no extra options - the code returns a mapcube only
    test_output = mapcube_coalign_by_match_template(mc)
    assert(isinstance(test_output, map.MapCube))

    # Test returning with no clipping.  Output layers should have the same size
    # as the original input layer.
    test_mc = mapcube_coalign_by_match_template(mc, clip=False)
    assert(test_mc[0].data.shape == testmap.data.shape)
    assert(test_mc[1].data.shape == testmap.data.shape)

