"""Test cases for SOHO Map subclasses.
This particular test file pertains to LASCOMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

from sunpy.map.sources.soho import LASCOMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "lasco_c2_25299383_s.fts"))
lasco = Map(fitspath)

# LASCO Tests
def test_fitstoEIT():
    """Tests the creation of LASCOMap using FITS."""
    assert isinstance(lasco, LASCOMap)

def test_is_datasource_for():
    """Test the is_datasource_for method of LASCOMap.
    Note that header data to be provided as an argument
    can be a MapMeta object."""
    assert lasco.is_datasource_for(lasco.data, lasco.meta)

def test_measurement():
    """Tests the measurement property of the LASCOMap object."""
    assert lasco.measurement == "white-light"

def test_observatory():
    """Tests the observatory property of the LASCOMap object."""
    assert lasco.observatory == "SOHO"
