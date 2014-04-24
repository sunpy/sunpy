"""Test cases for SDO Map subclasses.
This particular test file pertains to AIAMap.
@Author: Pritish C. (VaticanCameos)
"""

import pytest

from sunpy.map.sources.sdo import AIAMap
from sunpy.map import Map
from sunpy import AIA_171_IMAGE as aiaimg

@pytest.fixture
def createAIAMap():
    """Creates an AIAMap as given in documentation examples, through AIA_171_IMAGE."""
    aia = Map(aiaimg)
    return aia

# AIA Tests
def test_AIAMap():
    """Tests the creation of AIAMap from AIA_171_IMAGE."""
    assert (isinstance(Map(aiaimg), AIAMap) == True)

def test_is_datasource_for(createAIAMap):
    """Tests the is_datasource_for method of AIAMap."""
    assert (createAIAMap.is_datasource_for(createAIAMap.data, createAIAMap.meta) == True)

def test_observatory(createAIAMap):
    """Tests the observatory property of the AIAMap object."""
    assert (createAIAMap.observatory == "SDO")

def test_measurement(createAIAMap):
    """Tests the measurement property of the AIAMap object."""
    assert (createAIAMap.measurement == 171)
