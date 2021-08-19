"""Test cases for SJIMap.
This particular test file pertains to SJIMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

import pytest

import astropy.units as u

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.iris import SJIMap
from sunpy.util.exceptions import SunpyMetadataWarning, SunpyUserWarning


@pytest.fixture
def irismap():
    path = sunpy.data.test.rootdir
    fitspath = glob.glob(os.path.join(
        path, "iris_l2_20130801_074720_4040000014_SJI_1400_t000.fits"))
    with pytest.warns(SunpyUserWarning, match='This file contains more than 2 dimensions'):
        return Map(fitspath, silence_errors=True)


# IRIS Tests
def test_fitstoIRIS(irismap):
    """Tests the creation of SJIMap using FITS."""
    assert (isinstance(irismap, SJIMap))


def test_is_datasource_for(irismap):
    """Test the is_datasource_for method of SJIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert irismap.is_datasource_for(irismap.data, irismap.meta)


def test_observatory(irismap):
    """Tests the observatory property of SJIMap."""
    assert irismap.observatory == "IRIS"


def test_wcs(irismap):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        irismap.pixel_to_world(0*u.pix, 0*u.pix)
