"""
Test cases for SUVI Map subclass.
"""

import os
import glob

import pytest

from sunpy.map.sources.suvi import SUVIMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "test_suvi_image.fits"))


@pytest.fixture()
def suvi():
    """Creates an SUVIMap from a FITS file."""
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
