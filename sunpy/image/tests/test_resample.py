# Author: Tomas Meszaros <exo@tty.sk>

from sunpy import AIA_171_IMAGE
from sunpy import map
from sunpy.image.rescale import reshape_image_to_4d_superpixel


AIA_MAP = map.Map(AIA_171_IMAGE)

def resample_meta(dimensions, method, center, minusone):
	map_resampled = AIA_MAP.resample(dimensions)
	return map_resampled.shape

def resample_method(method):
	assert resample_meta((512, 512), method, False, False) == (512, 512)
	assert resample_meta((2056, 2056), method, False, False) == (2056, 2056)
	assert resample_meta((512, 512), method, False, True) == (512, 512)
	assert resample_meta((2056, 2056), method, False, True) == (2056, 2056)
	assert resample_meta((512, 512), method, True, False) == (512, 512)
	assert resample_meta((2056, 2056), method, True, False) == (2056, 2056)
	assert resample_meta((512, 512), method, True, True) == (512, 512)
	assert resample_meta((2056, 2056), method, True, True) == (2056, 2056)

def test_resample_neighbor():
	resample_method('neighbor')

def test_resample_nearest():
	resample_method('nearest')

def test_resample_linear():
	resample_method('linear')

def test_resample_spline():
	resample_method('spline')

def test_reshape():
	assert reshape_image_to_4d_superpixel(AIA_MAP.data, (512, 512)) != None
	assert reshape_image_to_4d_superpixel(AIA_MAP.data, (600, 512)) == None
	assert reshape_image_to_4d_superpixel(AIA_MAP.data, (512, 600)) == None

