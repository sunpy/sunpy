import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.suvi import SUVIMap


@pytest.fixture
def suvi_map():
    path = get_test_filepath("dr_suvi-l2-ci195_g16_s20190403T093200Z_e20190403T093600Z_v1-0-0_rebinned.header")
    return get_dummy_map_from_header(path)


def test_suvimap_creation(suvi_map):
    assert isinstance(suvi_map, SUVIMap)


def test_is_datasource_for(suvi_map):
    assert suvi_map.is_datasource_for(suvi_map.data, suvi_map.meta)


def test_observatory(suvi_map):
    assert suvi_map.observatory == "GOES-R"


def test_detector(suvi_map):
    assert suvi_map.detector == "SUVI"


def test_norm_clip(suvi_map):
    # Tests that the default normalizer has clipping disabled
    assert not suvi_map.plot_settings['norm'].clip


# SUVI provides observer coordinate information in an OBSGEO system, so this test
# needs remote data to access the latest IERS table to do a coordinate transformation from
# OBSGEO to heliographic Stonyhurst coordinates.
@pytest.mark.remote_data
def test_wcs(suvi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    suvi_map.pixel_to_world(0*u.pix, 0*u.pix)
