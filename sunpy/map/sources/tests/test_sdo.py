"""Test cases for SDO Map subclasses"""

from sunpy.map.sources.sdo import HMIMap
from sunpy.map.sources.sdo import AIAMap
from sunpy.map import Map
from sunpy.net.helioviewer import HelioviewerClient
import numpy as np
import pytest

@pytest.fixture
def createHMIMap():
    hv = HelioviewerClient()
    filepath = hv.download_jp2('2012/07/05 00:30:00', observatory='SDO', instrument='HMI', detector='HMI', measurement='continuum')
    hmi = Map(filepath)
    return hmi

def test_is_datasource_for(createHMIMap):
    assert (createHMIMap.is_datasource_for() == True)
