"""
Test cases for PUNCH subclass.
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.punch import PUNCHMap


@pytest.fixture
def punch_map():
    return get_dummy_map_from_header(get_test_filepath("punch.header"))


def test_punch_map(punch_map):
    """Tests the creation of PUNCHMap"""
    assert isinstance(punch_map, PUNCHMap)


def test_reference_date(punch_map):
    """ Tests the reference_date"""
    assert punch_map.reference_date.isot == "2025-03-12T03:26:00.000"


def test_date(punch_map):
    assert punch_map.date.isot == "2025-03-12T03:26:00.000"


def test_unit(punch_map):
    assert punch_map.unit == u.Unit("2.009e+07 W / (sr m2)")


def test_is_datasource_for(punch_map):
    """Tests the is_datasource_for method of PUNCHMap."""
    assert punch_map.is_datasource_for(punch_map.data, punch_map.meta)


def test_observatory(punch_map):
    """Tests the observatory property of the PUNCHMap object."""
    assert punch_map.observatory == "PUNCH"


def test_instrument(punch_map):
    """Tests the instrument property of the PUNCHMap object."""
    assert punch_map.instrument == "WFI+NFI Mosaic"


def test_wcs(punch_map):
    """Test WCS is valid and can transform from pixels to world coordinates"""
    punch_map.pixel_to_world(0*u.pix, 0*u.pix)
