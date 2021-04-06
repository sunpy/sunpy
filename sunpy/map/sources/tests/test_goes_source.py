import os
import glob

import pytest

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.goes import SXIMap

path = sunpy.data.test.rootdir


@pytest.fixture
def sxi():
    # PTHNA filtered image
    fitspath = glob.glob(os.path.join(path, "SXI_20170201_000400155_BA_13.FTS"))
    return Map(fitspath)


def test_fitstoSXI(sxi):
    """
    Tests the creation of SXIMap to fits.
    """
    assert isinstance(sxi, SXIMap)


def test_is_datasource_for(sxi):
    """
    Test the is_data_source_for method of SXIMap.
    """
    assert sxi.is_datasource_for(sxi.data, sxi.meta)


def test_measurement(sxi):
    """
    Tests the measurement property of SXIMap.
    """
    assert sxi.measurement == "PTHNA"


def test_observatory(sxi):
    """
    Tests the observatory property of SXIMap.
    """
    assert sxi.observatory == "GOES-13"


def test_norm_clip(sxi):
    """
    Tests that the default normalizer has clipping disabled.
    """
    assert not sxi.plot_settings['norm'].clip


def test_rotate(sxi):
    """
    Tests the rotate works (at least doesn't error).
    """
    rot_map = sxi.rotate()
    assert "pc1_1" in rot_map.meta
    assert "pc1_2" in rot_map.meta
    assert "pc2_1" in rot_map.meta
    assert "pc2_2" in rot_map.meta
