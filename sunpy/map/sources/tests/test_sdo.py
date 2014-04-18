"""Test cases for SDO Map subclasses"""

import sunpy
from sunpy.map import HMIMap
from sunpy.map import AIAMap
from sunpy.map import Map
from sunpy.net.helioviewer import HelioviewerClient
import numpy as np


@pytest.fixture
def createHMIMap():
    hv = HelioviewerClient()
    filepath = hv.download_jp2('2012/07/05 00:30:00', observatory='SDO', instrument='HMI', detector='HMI', measurement='continuum')
    hmi = Map(filepath)
    return hmi

def test_is_datasource_for(createHMIMap):
    assert (
