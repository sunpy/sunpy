"""Test cases for SOHO Map subclasses.
This particular test file pertains to LASCOMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

import pytest

from sunpy.map.sources.soho import LASCOMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "lasco_c2_25299383_s.fts"))

@pytest.fixture
def createLASCO():
    """Creates a LASCOMap from a FITS file."""
    return Map(fitspath)

# LASCO Tests
def test_fitstoEIT(createLASCO):
    """Tests the creation of LASCOMap using FITS."""
    assert (isinstance(createLASCO, LASCOMap) == True)

def test_is_datasource_for(createLASCO):
    """Test the is_datasource_for method of LASCOMap.
    Note that header data to be provided as an argument
    can be a MapMeta object."""
    assert (createLASCO.is_datasource_for(createLASCO.data, createLASCO.meta) == True)

def test_measurement(createLASCO):
    """Tests the measurement property of the LASCOMap object."""
    assert (createLASCO.measurement == "white-light")

def test_observatory(createLASCO):
    """Tests the observatory property of the LASCOMap object."""
    assert(createLASCO.observatory == "SOHO")
