import os

import pytest
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u

import sunpy.data.test
import sunpy.map
from sunpy.physics import limb_darkening


@pytest.fixture
def aia171_test_map():
    testpath = sunpy.data.test.rootdir
    return sunpy.map.Map(os.path.join(testpath, 'aia_171_level1.fits'))

def test_limb_filter(aia171_test_map):
    with pytest.raises(ValueError):
        limb_darkening._get_limb_dark(aia171_test_map)
