"""Test cases for SDO Map subclasses.
This particular test file pertains to AIAMap.
@Author: Pritish C. (VaticanCameos)
"""
import os
import glob

import pytest

from sunpy.map.sources.sdo import AIAMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
jp2path = glob.glob(os.path.join(path, "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))
aiaimg = glob.glob(os.path.join(path, "aia_171_level1.fits"))

from sunpy.tests.helpers import SKIP_GLYMUR

if SKIP_GLYMUR:
    params = [aiaimg]
else:
    params = [aiaimg, jp2path]

# The fixture is parameterized with aiaimg and jp2path.
@pytest.fixture(scope="module", params=params)
def createAIAMap(request):
    """Creates an AIAMap as given in documentation examples, through AIA_171_IMAGE
    or through the use of the JP2 file."""
    aiaobj = Map(request.param)
    return aiaobj

# AIA Tests
def test_AIAMap(createAIAMap):
    """Tests the creation of AIAMap from AIA_171_IMAGE or through
    use of the JP2 file."""
    assert isinstance(createAIAMap, AIAMap)

def test_is_datasource_for(createAIAMap):
    """Tests the is_datasource_for method of AIAMap."""
    assert createAIAMap.is_datasource_for(createAIAMap.data, createAIAMap.meta)

def test_observatory(createAIAMap):
    """Tests the observatory property of the AIAMap object."""
    assert createAIAMap.observatory == "SDO"

def test_measurement(createAIAMap):
    """Tests the measurement property of the AIAMap object."""
    assert createAIAMap.measurement.value in [171, 193]
    # aiaimg has 171, jp2path has 193.
