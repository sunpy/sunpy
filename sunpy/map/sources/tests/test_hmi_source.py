import os
import glob

import pytest

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.sdo import HMIMap


@pytest.fixture
def hmimap():
    path = sunpy.data.test.rootdir
    fitspath = glob.glob(os.path.join(path, "resampled_hmi.fits"))
    return Map(fitspath)


def test_fitstoHMI(hmimap):
    """
    Tests the creation of HMIMap using FITS.
    """
    assert isinstance(hmimap, HMIMap)


def test_is_datasource_for(hmimap):
    """
    Test the is_datasource_for method of HMIMap.
    """
    assert hmimap.is_datasource_for(hmimap.data, hmimap.meta)


def test_observatory(hmimap):
    """
    Tests the observatory property of HMIMap.
    """
    assert hmimap.observatory == "SDO"


def test_measurement(hmimap):
    """
    Tests the measurement property of HMIMap.
    """
    assert hmimap.measurement == "continuum"
