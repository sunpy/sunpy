from sunpy import AIA_171_IMAGE
from sunpy import map
from sunpy.image.normalization import multiscale_gaussian


def test_gaussian():
    data = map.Map(AIA_171_IMAGE).data
    new_image = multiscale_gaussian(data)
    assert new_image.shape == data.shape