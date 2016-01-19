"""Test cases for SDO Map subclasses.
This particular test file pertains to HMIMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

from sunpy.map.sources.sdo import HMIMap
from sunpy.map import Map
import sunpy.data.test
#from sunpy.net import HelioviewerClient

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "resampled_hmi.fits"))
hmi = Map(fitspath)

# HMI Tests
def test_fitstoHMI():
    """Tests the creation of HMIMap using FITS."""
    assert isinstance(hmi, HMIMap)

def test_is_datasource_for():
    """Test the is_datasource_for method of HMIMap.
    Note that header data to be provided as an argument
    can be a MapMeta object, which in this case is
    hmi.meta."""
    assert hmi.is_datasource_for(hmi.data, hmi.meta)

def test_observatory():
    """Tests the observatory property of the HMIMap object."""
    assert hmi.observatory == "SDO"

def test_measurement():
    """Tests the measurement property of the HMIMap object."""
    assert hmi.measurement == "continuum"
