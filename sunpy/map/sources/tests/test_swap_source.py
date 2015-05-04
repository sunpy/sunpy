"""Test cases for PROBA2 Map subclasses.
This particular test file pertains to SWAPMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

import pytest

from sunpy.map.sources.proba2 import SWAPMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
fitslist = glob.glob(os.path.join(path, "SWAP", "*"))

@pytest.fixture(scope="module", params=fitslist)
def createSWAP(request):
    """Creates an SWAPMap from a FITS file."""
    return Map(request.param)

# SWAP Tests
def test_fitstoSWAP(createSWAP):
    """Tests the creation of SWAPMap using FITS."""
    assert isinstance(createSWAP, SWAPMap)

def test_is_datasource_for(createSWAP):
    """Test the is_datasource_for method of SWAPMap.
    Note that header data to be provided as an argument
    can be a MapMeta object."""
    assert createSWAP.is_datasource_for(createSWAP.data, createSWAP.meta)

def test_observatory(createSWAP):
    """Tests the observatory property of the SWAPMap object."""
    assert createSWAP.observatory == "PROBA2"

def test_measurement(createSWAP):
    """Tests the measurement property of the SWAPMap object."""
    assert createSWAP.measurement.value == 174
