"""
Test cases for STEREO EUVIMap subclass.
"""
import copy

import pytest

import astropy.units as u

from sunpy.coordinates import sun
from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.stereo import EUVIMap
from sunpy.sun import constants
from .helpers import _test_private_date_setters

__all__ = "Pritish C. (VaticanCameos)"


@pytest.fixture
def euvi_map():
    return get_dummy_map_from_header(get_test_filepath("euvi_20090615_000900_n4euA_s.header"))


def test_fitstoEUVI(euvi_map):
    """Tests the creation of EUVIMap using FITS."""
    assert isinstance(euvi_map, EUVIMap)


def test_reference_date(euvi_map):
    assert euvi_map.reference_date.isot == "2009-06-15T00:09:08.009"


def test_date(euvi_map):
    assert euvi_map.date.isot == "2009-06-15T00:09:00.006"


def test_private_date_setters(euvi_map):
    _test_private_date_setters(euvi_map)


def test_is_datasource_for(euvi_map):
    """Test the is_datasource_for method of EUVIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert euvi_map.is_datasource_for(euvi_map.data, euvi_map.meta)


def test_measurement(euvi_map):
    """Tests the measurement property of the EUVIMap object."""
    assert euvi_map.measurement.value == 171


def test_observatory(euvi_map):
    """Tests the observatory property of the EUVIMap object."""
    assert euvi_map.observatory == "STEREO A"


def test_rsun_obs(euvi_map):
    """Tests the rsun_obs property"""
    assert euvi_map.rsun_obs.value == euvi_map.meta['rsun']


def test_rsun_missing(euvi_map):
    """Tests output if 'rsun' is missing"""
    euvi_no_rsun = euvi_map._new_instance(euvi_map.data, copy.deepcopy(euvi_map.meta))
    euvi_no_rsun.meta.pop('rsun', None)
    r = euvi_no_rsun.observer_coordinate.radius
    assert euvi_no_rsun.rsun_obs == sun._angular_radius(constants.radius, r)


def test_norm_clip(euvi_map):
    # Tests that the default normalizer has clipping disabled
    assert not euvi_map.plot_settings['norm'].clip


def test_wcs(euvi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    euvi_map.pixel_to_world(0*u.pix, 0*u.pix)
