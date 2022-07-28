
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.yohkoh import SXTMap

__author__ = "It's a me, Mario!"


@pytest.fixture
def sxt_map():
    return get_dummy_map_from_header(get_test_filepath("sxtf_20011113_005113_121.header"))


def test_sxt_map(sxt_map):
    assert isinstance(sxt_map, SXTMap)


def test_is_datasource_for_sxt(sxt_map):
    assert sxt_map.is_datasource_for(sxt_map.data, sxt_map.meta)


def test_measurement_sxt(sxt_map):
    assert sxt_map.measurement == "Al01"


def test_observatory_sxt(sxt_map):
    assert sxt_map.observatory == "Yohkoh"


def test_rsun_obs_sxt(sxt_map):
    assert sxt_map.rsun_obs.value == sxt_map.meta['SOLAR_R']


def test_norm_clip_sxt(sxt_map):
    # Tests that the default normalizer has clipping disabled
    assert not sxt_map.plot_settings['norm'].clip


def test_wcs_sxt(sxt_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    sxt_map.pixel_to_world(0*u.pix, 0*u.pix)
