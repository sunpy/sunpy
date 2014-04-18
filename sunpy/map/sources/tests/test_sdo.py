"""Test cases for SDO Map subclasses"""

from sunpy.map.sources.sdo import HMIMap
from sunpy.map.sources.sdo import AIAMap
from sunpy.map import Map
from sunpy.io import jp2 as jp
from sunpy.net.helioviewer import HelioviewerClient
import pytest

@pytest.fixture
def createHMI():
    hv = HelioviewerClient()
    filepath = hv.download_jp2('2012/07/05 00:30:00', observatory='SDO', instrument='HMI', detector='HMI', measurement='continuum')
    return filepath

def test_is_datasource_for(createHMI):
    hmi = Map(createHMI)
    header = jp.get_header(createHMI)
    assert (hmi.is_datasource_for(hmi.data, header) == True)

def test_observatory(createHMI):
    hmi = Map(createHMI)
    assert(hmi.observatory == "SDO")

def test_measurement(createHMI):
    hmi = Map(createHMI)
    assert(hmi.measurement == "continuum")
    
