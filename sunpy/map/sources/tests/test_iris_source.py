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
irismap = Map(fitspath)


# IRIS Tests
def test_fitstoIRIS():
    """Tests the creation of SJIMap using FITS."""
    assert (isinstance(irismap, SJIMap))


def test_is_datasource_for():
    """Test the is_datasource_for method of SJIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert irismap.is_datasource_for(irismap.data, irismap.meta)


def test_observatory():
    """Tests the observatory property of SJIMap."""
    assert irismap.observatory == "IRIS"
