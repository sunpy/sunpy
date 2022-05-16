"""
Test cases for PROBA2 SWAPMap subclass.
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.proba2 import SWAPMap

header_files = [
    'SWAP/resampled0_swap.header',
    'SWAP/resampled1_swap.header',
    'SWAP/resampled2_swap.header',
    'SWAP/resampled3_swap.header',
]

__author__ = 'Pritish C. (VaticanCameos)'


@pytest.fixture(scope="module", params=header_files)
def swap_map(request):
    return get_dummy_map_from_header(get_test_filepath(request.param))


def test_fitstoSWAP(swap_map):
    """Tests the creation of SWAPMap using FITS."""
    assert isinstance(swap_map, SWAPMap)


def test_is_datasource_for(swap_map):
    """Test the is_datasource_for method of SWAPMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert swap_map.is_datasource_for(swap_map.data, swap_map.meta)


def test_observatory(swap_map):
    """Tests the observatory property of the SWAPMap object."""
    assert swap_map.observatory == "PROBA2"


def test_measurement(swap_map):
    """Tests the measurement property of the SWAPMap object."""
    assert swap_map.measurement.value == 174


def test_wcs(swap_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    swap_map.pixel_to_world(0*u.pix, 0*u.pix)
