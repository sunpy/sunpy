"""Test cases for SDO Map subclasses.
This particular test file pertains to HMIMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

import pytest

from sunpy.map.sources.sdo import HMIMap
from sunpy.map import Map
import sunpy.data.test
#from sunpy.net import HelioviewerClient

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "resampled_hmi.fits"))
hmi = Map(fitspath)

# This fixture is no longer used on Stuart's instructions. Downloading is expensive.
@pytest.mark.online
@pytest.fixture
def createHMI():
    """Downloads a HMIMap jp2 object through the use of HelioviewerClient."""
    hv = HelioviewerClient()
    filepath = hv.download_jp2('2012/07/05 00:30:00', observatory='SDO', instrument='HMI', detector='HMI', measurement='continuum')
    return filepath

# HMI Tests
def test_fitstoHMI():
    """Tests the creation of HMIMap using FITS."""
    assert (isinstance(hmi, HMIMap) == True)
    
def test_is_datasource_for():
    """Test the is_datasource_for method of HMIMap.
    Note that header data to be provided as an argument
    can be a MapMeta object, which in this case is
    hmi.meta."""
    assert (hmi.is_datasource_for(hmi.data, hmi.meta) == True)

def test_observatory():
    """Tests the observatory property of the HMIMap object."""
    assert(hmi.observatory == "SDO")

def test_measurement():
    """Tests the measurement property of the HMIMap object."""
    assert (hmi.measurement == "continuum")

