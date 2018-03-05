import os
import glob

from sunpy.map.sources.trace import TRACEMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "tsi20010130_025823_a2.fits"))
trace = Map(fitspath)

# TRACE Tests
def test_fitstoTRACE():
    """Tests the creation of TRACEMap using FITS."""
    assert isinstance(trace, TRACEMap)

def test_is_datasource_for():
    """Test the is_datasource_for method of TRACEMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert trace.is_datasource_for(trace.data, trace.meta)

def test_measurement():
    """Tests the measurement property of the TRACEMap object."""
    assert int(trace.measurement) == 171

def test_observatory():
    """Tests the observatory property of the TRACEMap object."""
    assert trace.observatory == "TRACE"