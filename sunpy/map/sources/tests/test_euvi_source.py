"""Test cases for STEREO Map subclasses.
This particular test file pertains to EUVIMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

import pytest

import sunpy.data.test
from sunpy.coordinates import sun
from sunpy.map import Map
from sunpy.map.sources.stereo import EUVIMap
from sunpy.sun import constants
from sunpy.util.exceptions import SunpyUserWarning

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "euvi_20090615_000900_n4euA_s.fts"))
euvi = Map(fitspath)


# EUVI Tests
def test_fitstoEUVI():
    """Tests the creation of EUVIMap using FITS."""
    assert isinstance(euvi, EUVIMap)


def test_is_datasource_for():
    """Test the is_datasource_for method of EUVIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert euvi.is_datasource_for(euvi.data, euvi.meta)


def test_measurement():
    """Tests the measurement property of the EUVIMap object."""
    assert euvi.measurement.value == 171


def test_observatory():
    """Tests the observatory property of the EUVIMap object."""
    assert euvi.observatory == "STEREO A"


def test_rsun_obs():
    """Tests the rsun_obs property"""
    assert euvi.rsun_obs.value == euvi.meta['rsun']


def test_rsun_missing():
    """Tests output if 'rsun' is missing"""
    euvi_no_rsun = Map(fitspath)
    euvi_no_rsun.meta['rsun'] = None
    r = euvi_no_rsun.observer_coordinate.radius
    with pytest.warns(SunpyUserWarning, match='Missing metadata for solar angular radius'):
        assert euvi_no_rsun.rsun_obs == sun._angular_radius(constants.radius, r)


def test_norm_clip():
    # Tests that the default normalizer has clipping disabled
    assert not euvi.plot_settings['norm'].clip
