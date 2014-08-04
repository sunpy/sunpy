from sunpy import AIA_171_IMAGE
from sunpy import map
from sunpy.image.normalization import multiscale_gaussian


def test_gaussian():
    data = map.Map(AIA_171_IMAGE).data
    new_image = multiscale_gaussian(data)
    assert new_image.shape == data.shape


def test_gaussian_weights():
    data = map.Map(AIA_171_IMAGE).data
    new_image = multiscale_gaussian(data, weights=[1, 2, 3, 4, 5, 6])
    assert new_image.shape == data.shape