"""Test cases for STEREO Map subclasses.
This particular test file pertains to EUVIMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob
import sunpy.data.test

from sunpy.map import Map
from sunpy.map.sources.stereo import EUVIMap
from sunpy.map.sources.stereo import CORMap

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "euvi_20090615_000900_n4euA_s.fts"))
euvi = Map(fitspath)
CORFITSPATH = glob.glob(os.path.join(path, "cor1_20090615_000500_s4c1A.fts"))
CORMAP = Map(CORFITSPATH)

# EUVI Tests
def test_fitstoEIT():
    """Tests the creation of EUVIMap using FITS."""
    assert isinstance(euvi, EUVIMap)

def test_is_datasource_for():
    """Test the is_datasource_for method of EUVIMap.
    Note that header data to be provided as an argument
    can be a MapMeta object."""
    assert euvi.is_datasource_for(euvi.data, euvi.meta)

def test_measurement():
    """Tests the measurement property of the EUVIMap object."""
    assert euvi.measurement == 171

def test_observatory():
    """Tests the observatory property of the EUVIMap object."""
    assert euvi.observatory == "STEREO_A"

def test_rsun_arcseconds():
    """ Tests the rsun_arcseconds property. """
    assert euvi.rsun_arcseconds == euvi.meta['rsun']

#CORMap Tests
def test_CORMap_measurement():
    """Tests the measurement property of the CORMap object."""
    assert CORMAP.measurement == "white-light"
