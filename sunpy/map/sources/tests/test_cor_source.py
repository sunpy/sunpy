"""Test cases for STEREO Map subclasses.
This particular test file pertains to CORMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.stereo import CORMap

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "cor1_20090615_000500_s4c1A.fts"))
cor = Map(fitspath)

# COR Tests


def test_fitstoCOR():
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


def test_norm_clip():
    # Tests that the default normalizer has clipping disabled
    assert not cor.plot_settings['norm'].clip
