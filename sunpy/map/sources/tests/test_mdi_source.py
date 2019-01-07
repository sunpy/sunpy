"""Test cases for SOHO Map subclasses.
This particular test file pertains to MDIMap.
@Author: Pritish C. (VaticanCameos)
"""

<<<<<<< HEAD
from pathlib import Path
=======
import os
import pathlib
>>>>>>> pathlib added astropy files not touched yet
import glob

from sunpy.map.sources.soho import MDIMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
<<<<<<< HEAD
fitspath = glob.glob(str(Path.home().joinpath(path, "mdi_fd_Ic_6h_01d.5871.0000_s.fits")))
=======
fitspath = glob.glob(str(pathlib.Path.home().joinpath(path, "mdi_fd_Ic_6h_01d.5871.0000_s.fits")))
>>>>>>> pathlib added astropy files not touched yet
mdi = Map(fitspath)

# MDI Tests
def test_fitstoMDI():
    """Tests the creation of MDIMap using FITS."""
    assert isinstance(mdi, MDIMap)

def test_is_datasource_for():
    """Test the is_datasource_for method of MDIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert mdi.is_datasource_for(mdi.data, mdi.meta)

def test_observatory():
    """Tests the observatory property of the MDIMap object."""
    assert mdi.observatory == "SOHO"

def test_measurement():
    """Tests the measurement property of the MDIMap object."""
    assert mdi.measurement == "continuum"
