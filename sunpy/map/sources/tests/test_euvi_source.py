import os
import glob

import pytest

import sunpy.data.test
from sunpy.coordinates import sun
from sunpy.map import Map
from sunpy.map.sources.stereo import EUVIMap
from sunpy.sun import constants
from sunpy.util.exceptions import SunpyUserWarning


@pytest.fixture
def euvimap():
    path = sunpy.data.test.rootdir
    fitspath = glob.glob(os.path.join(path, "euvi_20090615_000900_n4euA_s.fts"))
    return Map(fitspath)


def test_fitstoEUVI(euvimap):
    """
    Tests the creation of EUVIMap using FITS.
    """
    assert isinstance(euvimap, EUVIMap)


def test_is_datasource_for(euvimap):
    """
    Test the is_datasource_for method of EUVIMap.
    """
    assert euvimap.is_datasource_for(euvimap.data, euvimap.meta)


def test_measurement(euvimap):
    """
    Tests the measurement property of EUVIMap.
    """
    assert euvimap.measurement.value == 171


def test_observatory(euvimap):
    """
    Tests the observatory property of EUVIMap.
    """
    assert euvimap.observatory == "STEREO A"


def test_rsun_obs(euvimap):
    """
    Tests the rsun_obs property.
    """
    assert euvimap.rsun_obs.value == euvimap.meta['rsun']


def test_rsun_missing(euvimap):
    """
    Tests output if 'rsun' is missing.
    """
    euvimap.meta['rsun'] = None
    r = euvimap.observer_coordinate.radius
    with pytest.warns(SunpyUserWarning, match='Missing metadata for solar angular radius'):
        assert euvimap.rsun_obs == sun._angular_radius(constants.radius, r)


def test_norm_clip(euvimap):
    """
    Tests that the default normalizer has clipping disabled.
    """
    assert not euvimap.plot_settings['norm'].clip
