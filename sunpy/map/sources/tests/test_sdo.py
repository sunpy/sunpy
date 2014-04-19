"""Test cases for SDO Map subclasses"""

from sunpy.map.sources.sdo import HMIMap
from sunpy.map.sources.sdo import AIAMap
from sunpy.map import Map
from sunpy.net import HelioviewerClient
from sunpy import AIA_171_IMAGE as aiaimg
import pytest

@pytest.mark.online
@pytest.fixture
def createHMI():
    """Downloads a HMIMap jp2 object through the use of HelioviewerClient."""
    hv = HelioviewerClient()
    filepath = hv.download_jp2('2012/07/05 00:30:00', observatory='SDO', instrument='HMI', detector='HMI', measurement='continuum')
    return filepath

@pytest.fixture
def createAIAMap():
    """Creates an AIAMap as given in documentation examples, through AIA_171_IMAGE."""
    aia = Map(aiaimg)
    return aia

# HMI Tests
def test_is_datasource_for(createHMI):
    hmi = Map(createHMI)
    header = dict(hmi.meta)
    assert (hmi.is_datasource_for(hmi.data, header) == True)

def test_observatory(createHMI):
    hmi = Map(createHMI)
    assert(hmi.observatory == "SDO")

def test_measurement(createHMI):
    hmi = Map(createHMI)
    assert (hmi.measurement == "continuum")

# AIA Tests
def test_is_datasource_for(createAIAMap):
    header = dict(createAIAMap.meta)
    assert (createAIAMap.is_datasource_for(createAIAMap.data, header) == True)

def test_observatory(createAIAMap):
    assert (createAIAMap.observatory == "SDO")

def test_measurement(createAIAMap):
    assert (createAIAMap.measurement == 171)
