import os

import pytest

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.trace import TRACEMap


@pytest.fixture
def createTRACE():
    """
    Creates a TRACEMap from a FITS file.
    """
    path = sunpy.data.test.rootdir
    fitspath = os.path.join(path, "tsi20010130_025823_a2.fits")
    return Map(fitspath)


def test_fitstoTRACE(createTRACE):
    """
    Tests the creation of TRACEMap using FITS.
    """
    assert isinstance(createTRACE, TRACEMap)


def test_is_datasource_for(createTRACE):
    """
    Test the is_datasource_for method of TRACEMap.
    """
    assert createTRACE.is_datasource_for(createTRACE.data, createTRACE.meta)


def test_measurement(createTRACE):
    """
    Tests the measurement property of TRACEMap.
    """
    assert int(createTRACE.measurement) == 171


def test_observatory(createTRACE):
    """
    Tests the observatory property of TRACEMap.
    """
    assert createTRACE.observatory == "TRACE"


def test_norm_clip(createTRACE):
    """
    Tests that the default normalizer has clipping disabled.
    """
    assert not createTRACE.plot_settings['norm'].clip
