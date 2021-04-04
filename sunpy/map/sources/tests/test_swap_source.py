import os
import glob

import pytest

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.proba2 import SWAPMap

path = sunpy.data.test.rootdir
fitslist = glob.glob(os.path.join(path, "SWAP", "*"))


@pytest.fixture(params=fitslist)
def createSWAP(request):
    """
    Creates an SWAPMap from a FITS file.
    """
    return Map(request.param)


def test_fitstoSWAP(createSWAP):
    """
    Tests the creation of SWAPMap using FITS.
    """
    assert isinstance(createSWAP, SWAPMap)


def test_is_datasource_for(createSWAP):
    """
    Test the is_datasource_for method of SWAPMap.
    """
    assert createSWAP.is_datasource_for(createSWAP.data, createSWAP.meta)


def test_observatory(createSWAP):
    """
    Tests the observatory property of SWAPMap.
    """
    assert createSWAP.observatory == "PROBA2"


def test_measurement(createSWAP):
    """
    Tests the measurement property of SWAPMap.
    """
    assert createSWAP.measurement.value == 174
