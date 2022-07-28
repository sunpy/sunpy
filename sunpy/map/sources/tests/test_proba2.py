import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.proba2 import SWAPMap

__author__ = 'Pritish C. (VaticanCameos)'
swap_header_files = [
    'SWAP/resampled0_swap.header',
    'SWAP/resampled1_swap.header',
    'SWAP/resampled2_swap.header',
    'SWAP/resampled3_swap.header',
]


@pytest.fixture(scope="module", params=swap_header_files)
def swap_map(request):
    return get_dummy_map_from_header(get_test_filepath(request.param))


def test_swap_map(swap_map):
    assert isinstance(swap_map, SWAPMap)


def test_is_datasource_for_swap(swap_map):
    assert swap_map.is_datasource_for(swap_map.data, swap_map.meta)


def test_observatory_swap(swap_map):
    assert swap_map.observatory == "PROBA2"


def test_measurement_swap(swap_map):
    assert swap_map.measurement.value == 174


def test_wcs_swap(swap_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    swap_map.pixel_to_world(0*u.pix, 0*u.pix)
