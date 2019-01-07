"""Test cases for STEREO Map subclasses.
This particular test file pertains to CORMap.
@Author: Pritish C. (VaticanCameos)
"""

<<<<<<< HEAD
from pathlib import Path
=======
import os
import pathlib
>>>>>>> pathlib added astropy files not touched yet
import glob

from sunpy.map.sources.stereo import CORMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
<<<<<<< HEAD
fitspath = glob.glob(str(Path.home().joinpath(path, "cor1_20090615_000500_s4c1A.fts")))
=======
fitspath = glob.glob(str(pathlib.Path.home().joinpath(path, "cor1_20090615_000500_s4c1A.fts")))
>>>>>>> pathlib added astropy files not touched yet
cor = Map(fitspath)

# COR Tests
def test_fitstoEIT():
    """Tests the creation of CORMap using FITS."""
    assert isinstance(cor, CORMap)

def test_is_datasource_for():
    """Test the is_datasource_for method of CORMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert cor.is_datasource_for(cor.data, cor.meta)

def test_measurement():
    """Tests the measurement property of the CORMap object."""
    assert cor.measurement == "white-light"

def test_observatory():
    """Tests the observatory property of the CORMap object."""
    assert cor.observatory == "STEREO A"
