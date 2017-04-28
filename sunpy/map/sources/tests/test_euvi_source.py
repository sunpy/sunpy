"""Test cases for STEREO Map subclasses.
This particular test file pertains to EUVIMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

from sunpy.map.sources.stereo import EUVIMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "euvi_20090615_000900_n4euA_s.fts"))
euvi = Map(fitspath)

# EUVI Tests
def test_fitstoEIT():
    """Tests the creation of EUVIMap using FITS."""
    assert isinstance(euvi, EUVIMap)

def test_is_datasource_for():
    """Test the is_datasource_for method of EUVIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert euvi.is_datasource_for(euvi.data, euvi.meta)

def test_measurement():
    """Tests the measurement property of the EUVIMap object."""
    assert euvi.measurement.value == 171

def test_observatory():
    """Tests the observatory property of the EUVIMap object."""
    assert euvi.observatory == "STEREO A"
