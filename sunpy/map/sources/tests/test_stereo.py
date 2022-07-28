import copy

import pytest

import astropy.units as u

from sunpy.coordinates import sun
from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.stereo import CORMap, EUVIMap, HIMap
from sunpy.sun import constants


@pytest.fixture
def cor_map():
    header_file = get_test_filepath("cor1_20090615_000500_s4c1A.header")
    return get_dummy_map_from_header(header_file)


@pytest.fixture
def hi_map():
    return get_dummy_map_from_header(get_test_filepath('hi_20110910_114721_s7h2A.header'))


@pytest.fixture
def euvi_map():
    return get_dummy_map_from_header(get_test_filepath("euvi_20090615_000900_n4euA_s.header"))


def test_hi_map(hi_map):
    assert isinstance(hi_map, HIMap)


def test_is_datasource_for_hi(hi_map):
    assert hi_map.is_datasource_for(hi_map.data, hi_map.meta)


def test_measurement_hi(hi_map):
    assert hi_map.measurement == "white-light"


def test_observatory_hi(hi_map):
    assert hi_map.observatory == "STEREO A"


def test_norm_clip_hi(hi_map):
    # Tests that the default normalizer has clipping disabled
    assert not hi_map.plot_settings['norm'].clip


def test_wcs_hi(hi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    hi_map.pixel_to_world(0*u.pix, 0*u.pix)


def test_euvi_map(euvi_map):
    assert isinstance(euvi_map, EUVIMap)


def test_is_datasource_for_euvi(euvi_map):
    assert euvi_map.is_datasource_for(euvi_map.data, euvi_map.meta)


def test_measurement_euvi(euvi_map):
    assert euvi_map.measurement.value == 171


def test_observatory_euvi(euvi_map):
    assert euvi_map.observatory == "STEREO A"


def test_rsun_obs_euvi(euvi_map):
    assert euvi_map.rsun_obs.value == euvi_map.meta['rsun']


def test_rsun_missing_euvi(euvi_map):
    euvi_no_rsun = euvi_map._new_instance(euvi_map.data, copy.deepcopy(euvi_map.meta))
    euvi_no_rsun.meta.pop('rsun', None)
    r = euvi_no_rsun.observer_coordinate.radius
    assert euvi_no_rsun.rsun_obs == sun._angular_radius(constants.radius, r)


def test_norm_clip_euvi(euvi_map):
    # Tests that the default normalizer has clipping disabled
    assert not euvi_map.plot_settings['norm'].clip


def test_wcs_euvi(euvi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    euvi_map.pixel_to_world(0*u.pix, 0*u.pix)


def test_cor_map(cor_map):
    assert isinstance(cor_map, CORMap)


def test_is_datasource_for_cor(cor_map):
    assert cor_map.is_datasource_for(cor_map.data, cor_map.meta)


def test_measurement_cor(cor_map):
    assert cor_map.measurement == "white-light"


def test_observatory_cor(cor_map):
    assert cor_map.observatory == "STEREO A"


def test_norm_clip_cor(cor_map):
    # Tests that the default normalizer has clipping disabled
    assert not cor_map.plot_settings['norm'].clip


def test_wcs_cor(cor_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    cor_map.pixel_to_world(0*u.pix, 0*u.pix)
