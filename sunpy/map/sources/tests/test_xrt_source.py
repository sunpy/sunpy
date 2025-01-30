"""
Test cases for HINODE XRTMap subclass.
"""
import copy

import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.mapbase import SpatialPair
from sunpy.map.sources.hinode import XRTMap
from .helpers import _test_private_date_setters

__author__ = 'Pritish C. (VaticanCameos)'


@pytest.fixture
def xrt_map():
    return get_dummy_map_from_header(get_test_filepath("HinodeXRT.header"))


def test_fits_to_xrt(xrt_map):
    """Tests the creation of XRTMap using FITS."""
    assert isinstance(xrt_map, XRTMap)


def test_xrt_coordinate_system(xrt_map):
    assert xrt_map.coordinate_system ==  SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')


def test_is_datasource_for(xrt_map):
    """Test the is_datasource_for method of XRTMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert xrt_map.is_datasource_for(xrt_map.data, xrt_map.meta)


def test_reference_date(xrt_map):
    assert xrt_map.reference_date.isot == "2006-11-11T00:00:19.141"


def test_date(xrt_map):
    assert xrt_map.date.isot == "2006-11-11T00:00:19.141"


def test_private_date_setters(xrt_map):
    _test_private_date_setters(xrt_map)


def test_observatory(xrt_map):
    """Tests the observatory property of the XRTMap object."""
    assert xrt_map.observatory == "Hinode"


def test_measurement(xrt_map):
    """Tests the measurement property of the XRTMap object."""
    assert xrt_map.measurement == 'Be thin-Open'


def test_unit(xrt_map):
    """Tests the unit property of the XRTMap object."""
    assert xrt_map.unit == u.DN / u.second

def test_unit_no_renormalize(xrt_map):
    """
    Tests that the unit defaults ot DN if the history key does not say it has
    been normalized
    """
    new_xrt_map = copy.deepcopy(xrt_map)
    new_xrt_map.meta.pop('history')
    assert new_xrt_map.unit == u.DN

def test_level_number(xrt_map):
    assert xrt_map.processing_level == 1


def test_heliographic_longitude(xrt_map):
    assert u.allclose(xrt_map.heliographic_longitude, 0 * u.deg)


def test_heliographic_latitude(xrt_map):
    assert u.allclose(xrt_map.heliographic_latitude, 3.33047459 * u.deg)


def test_wheel_measurements(xrt_map):
    """Tests the filter_wheel_measurements objects present
    in the XRTMap object."""
    assert (xrt_map.filter_wheel1_measurements ==
            ["Al_med", "Al_poly", "Be_med", "Be_thin", "C_poly", "Open"])
    assert (xrt_map.filter_wheel2_measurements ==
            ["Open", "Al_mesh", "Al_thick", "Be_thick", "Gband", "Ti_poly"])


def test_wcs(xrt_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    xrt_map.pixel_to_world(0*u.pix, 0*u.pix)
