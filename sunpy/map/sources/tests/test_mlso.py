import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.mlso import KCorMap


@pytest.fixture
def kcor_map():
    return get_dummy_map_from_header(get_test_filepath("20181209_180305_kcor_l1.5_rebinned.header"))


def test_kcormap_creation(kcor_map):
    assert isinstance(kcor_map, KCorMap)


def test_is_datasource_for(kcor_map):
    assert kcor_map.is_datasource_for(kcor_map.data, kcor_map.meta)


def test_measurement(kcor_map):
    assert kcor_map.measurement == 735 * u.nm


def test_observatory(kcor_map):
    assert kcor_map.observatory == "MLSO"


def test_norm_clip(kcor_map):
    # Tests that the default normalizer has clipping disabled
    assert not kcor_map.plot_settings['norm'].clip


def test_wcs(kcor_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    kcor_map.pixel_to_world(0*u.pix, 0*u.pix)
