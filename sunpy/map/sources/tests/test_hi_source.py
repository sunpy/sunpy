"""Test cases for STEREO Map subclasses.
This particular test file pertains to HIMap.
"""

import os
import glob

from sunpy.map.sources.stereo import HIMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path,"hi_20110910_114721_s7h2A.fts"))
hi = Map(fitspath)

def test_fitstoHI():
    """Tests the creation of HIMap to fits"""
    assert isinstance(hi, HIMap)
        
def test_is_datasource_for():
    """Test the is_data_source_for method of HIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert hi.is_datasource_for(hi.data, hi.meta)

def test_measurement():
    """Tests the measurement property of the HIMap object."""
    assert hi.measurement == "white-light"

def test_observatory():
    """Tests the observatory property of the HIMap object."""
    assert hi.observatory == "STEREO A"

