"""Test cases for SJIMap.
This particular test file pertains to SJIMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

from sunpy.map import Map
from sunpy.map.sources.iris import SJIMap
from sunpy.map.mapbase import GenericMap
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "iris_l2_20130801_074720_4040000014_SJI_1400_t000.fits"))
irislist = Map(fitspath)

# IRIS Tests
def test_fitstoIRIS():
    """Tests the creation of SJIMap using FITS."""
    for amap in irislist:
        assert (isinstance(amap, (SJIMap, GenericMap)))

def test_is_datasource_for():
    """Test the is_datasource_for method of SJIMap.
    Note that header data to be provided as an argument
    can be a MapMeta object."""
    for amap in irislist:
        if isinstance(amap, SJIMap):
            assert amap.is_datasource_for(amap.data, amap.meta)

def test_observatory():
    """Tests the observatory property of SJIMap."""
    for amap in irislist:
        if isinstance(amap, SJIMap):
            assert amap.observatory == "IRIS"
