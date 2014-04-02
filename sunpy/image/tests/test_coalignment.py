# Author: Jack Ireland

import numpy as np
from numpy.testing import assert_allclose
from scipy.ndimage.interpolation import shift
from sunpy import AIA_171_IMAGE
from sunpy import map
from sunpy.image.coalignment import parabolic_turning_point, \
repair_nonfinite, default_fmap_function, _lower_clip, _upper_clip, \
calculate_clipping, get_correlation_shifts, find_best_match_location, \
match_template_to_layer, clip_edges, calculate_shift


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


def test_get_correlation_shifts():
	pass


def test_find_best_match_location():
	pass


def test_match_template_to_layer():
	pass

def test_lower_clip():
	print _lower_clip(clip_test_array)
	assert(_lower_clip(clip_test_array) == 2.0)

def test_upper_clip():
	print _upper_clip(clip_test_array)
	assert(_upper_clip(clip_test_array) == 1.0)

def test_calculate_clipping():
	answer = calculate_clipping(clip_test_array, clip_test_array)
	assert(answer == ([2.0, 1.0], [2.0, 1.0]))


def test_clip_edges():
	pass


def test_calculate_shift():
	pass


def test_default_fmap_function():
	assert(default_fmap_function([1,2,3]).dtype == np.float64(1).dtype)


def test_coalign():
	# take the AIA image and shift it
	m = map.Map(AIA_171_IMAGE)
	# Pixel displacements have the y-displacement as the first entry
	pixel_displacements = np.asarray([1.6, 10.1])
	known_displacements = {'x':np.asarray([0.0, pixel_displacements[1] * m.scale['x']]), 'y':np.asarray([0.0, pixel_displacements[0] * m.scale['y']])}
	# Apply the shift
	d1 = shift(m.data, pixel_displacements)
	m1 = map.Map((d1, m.meta))
	# Create the mapcube
	mc = map.Map([m, m1], cube=True)
	# Do the coalignment
	displacements = mc.coalign(displacements_only=True)
	# Assert
	assert_allclose(displacements['x'], known_displacements['x'], rtol=5e-2, atol=0)
	assert_allclose(displacements['y'], known_displacements['y'], rtol=5e-2, atol=0 )

