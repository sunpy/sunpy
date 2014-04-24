"""Test cases for SDO Map subclasses."""
"""This particular test file pertains to HMIMap."""

import pytest
from sunpy.map.sources.sdo import HMIMap
from sunpy.map import Map
#from sunpy.net import HelioviewerClient

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
    assert (isinstance(Map('resampled_hmi.fits'), HMIMap) == True)
    
def test_is_datasource_for():
    """Test the is_datasource_for method of HMIMap.
    Note that header data to be provided as an argument
    can be a MapMeta object, which in this case is
    hmi.meta."""
    hmi = Map('resampled_hmi.fits')
    assert (hmi.is_datasource_for(hmi.data, hmi.meta) == True)

def test_observatory():
    """Tests the observatory property of the HMIMap object."""
    hmi = Map('resampled_hmi.fits')
    assert(hmi.observatory == "SDO")

def test_measurement():
    """Tests the measurement property of the HMIMap object."""
    hmi = Map('resampled_hmi.fits')
    assert (hmi.measurement == "continuum")

