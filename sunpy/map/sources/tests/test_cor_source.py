import os
import glob

import pytest

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.stereo import CORMap


@pytest.fixture
def cormap():
    path = sunpy.data.test.rootdir
    fitspath = glob.glob(os.path.join(path, "cor1_20090615_000500_s4c1A.fts"))
    return Map(fitspath)


def test_fitstoCOR(cormap):
    """
    Tests the creation of CORMap using FITS.
    """
    assert isinstance(cormap, CORMap)


def test_is_datasource_for(cormap):
    """
    Test the is_datasource_for method of CORMap.
    """
    assert cormap.is_datasource_for(cormap.data, cormap.meta)


def test_measurement(cormap):
    """
    Tests the measurement property of CORMap.
    """
    assert cormap.measurement == "white-light"


def test_observatory(cormap):
    """
    Tests the observatory property of CORMap.
    """
    assert cormap.observatory == "STEREO A"


def test_norm_clip(cormap):
    """
    Tests that the default normalizer has clipping disabled.
    """
    assert not cormap.plot_settings['norm'].clip
