"""Test cases for SDO Map subclasses"""

import pytest
from sunpy.map.sources.sdo import HMIMap, AIAMap
from sunpy.map import Map
from sunpy.net import HelioviewerClient
from sunpy import AIA_171_IMAGE as aiaimg

@pytest.mark.online
@pytest.fixture
def createHMI():
    """Downloads a HMIMap jp2 object through the use of HelioviewerClient."""
    hv = HelioviewerClient()
    filepath = hv.download_jp2('2012/07/05 00:30:00', observatory='SDO', instrument='HMI', detector='HMI', 
    measurement='continuum')
    return filepath

@pytest.fixture
def createAIAMap():
    """Creates an AIAMap as given in documentation examples, through AIA_171_IMAGE."""
    return Map(aiaimg)

# HMI Tests
def test_is_datasource_for(createHMI):
    """Test the is_datasource_for method of HMIMap.
    Note that header data to be provided as an argument
    can be a MapMeta object, which in this case is
    hmi.meta."""
    hmi = Map(createHMI)
    header = dict(hmi.meta)
    assert hmi.is_datasource_for(hmi.data, header)

def test_observatory(createHMI):
    """Tests the observatory property of the HMIMap object."""
    hmi = Map(createHMI)
    assert hmi.observatory == "SDO"

def test_measurement(createHMI):
    """Tests the measurement property of the HMIMap object."""
    hmi = Map(createHMI)
    assert hmi.measurement == "continuum"
    
# AIA Tests
def test_is_datasource_for(createAIAMap):
    """Tests the is_datasource_for method of AIAMap."""
    header = dict(createAIAMap.meta)
    assert createAIAMap.is_datasource_for(createAIAMap.data, header)

def test_observatory(createAIAMap):
    """Tests the observatory property of the AIAMap object."""
    assert createAIAMap.observatory == "SDO"

def test_measurement(createAIAMap):
    """Tests the measurement property of the AIAMap object."""
    assert createAIAMap.measurement == 171
    
