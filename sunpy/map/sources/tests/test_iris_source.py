"""Test cases for IRISMap.
This particular test file pertains to IRISMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob
import numpy as np

import pytest

from sunpy.map import Map
from sunpy.map.sources.iris import IRISMap
from sunpy.map.mapbase import GenericMap
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "iris_l2_20130801_074720_4040000014_SJI_1400_t000.fits"))
irislist = Map(fitspath)

# IRIS Tests
def test_fitstoIRIS():
    """Tests the creation of IRISMap using FITS."""
    for amap in irislist:
        assert (isinstance(amap, (IRISMap, GenericMap)))

def test_is_datasource_for():
    """Test the is_datasource_for method of IRISMap.
    Note that header data to be provided as an argument
    can be a MapMeta object."""
    for amap in irislist:
        if isinstance(amap, IRISMap):
            assert amap.is_datasource_for(amap.data, amap.meta)

def test_observatory():
    """Tests the observatory property of IRISMap."""
    for amap in irislist:
        if isinstance(amap, IRISMap):
            assert amap.observatory == "IRIS"
