import os
import glob

import pytest

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.stereo import HIMap


@pytest.fixture
def himap():
    path = sunpy.data.test.rootdir
    fitspath = glob.glob(os.path.join(path, "hi_20110910_114721_s7h2A.fts"))
    return Map(fitspath)


def test_fitstoHI(himap):
    """
    Tests the creation of HIMap to fits.
    """
    assert isinstance(himap, HIMap)


def test_is_datasource_for(himap):
    """
    Test the is_data_source_for method of HIMap.
    """
    assert himap.is_datasource_for(himap.data, himap.meta)


def test_measurement(himap):
    """
    Tests the measurement property of HIMap.
    """
    assert himap.measurement == "white-light"


def test_observatory(himap):
    """
    Tests the observatory property of HIMap.
    """
    assert himap.observatory == "STEREO A"


def test_norm_clip(himap):
    """
    Tests that the default normalizer has clipping disabled.
    """
    assert not himap.plot_settings['norm'].clip
