"""
Test cases for SUVI Map subclass.
"""

import os
import glob

import pytest

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.suvi import SUVIMap


@pytest.fixture()
def suvi():
    """Creates an SUVIMap from a FITS file."""
    path = sunpy.data.test.rootdir
    fitspath = glob.glob(
        os.path.join(path, "dr_suvi-l2-ci195_g16_s20190403T093200Z_e20190403T093600Z_v1-0-0_rebinned.fits"))
    return Map(fitspath)


# SUVI Tests
def test_suvimap_creation(suvi):
    """Tests the creation of SUVIMap using FITS."""
    assert isinstance(suvi, SUVIMap)


def test_is_datasource_for(suvi):
    """Test the is_datasource_for method of SUVIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert suvi.is_datasource_for(suvi.data, suvi.meta)


def test_observatory(suvi):
    """Tests the observatory property of the SUVIMap object."""
    assert suvi.observatory == "GOES-R"


def test_detector(suvi):
    """Tests the detector property of the SUVIMap object."""
    assert suvi.detector == "SUVI"


def test_norm_clip(suvi):
    # Tests that the default normalizer has clipping disabled
    assert not suvi.plot_settings['norm'].clip
