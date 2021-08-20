"""
Test cases for HINODE XRTMap subclass.
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.hinode import XRTMap
from sunpy.util.exceptions import SunpyMetadataWarning

__author__ = 'Pritish C. (VaticanCameos)'


@pytest.fixture
def xrt():
    return get_dummy_map_from_header(get_test_filepath("HinodeXRT.header"))


# XRT Tests
def test_fitstoXRT(xrt):
    """Tests the creation of XRTMap using FITS."""
    assert isinstance(xrt, XRTMap)


def test_is_datasource_for(xrt):
    """Test the is_datasource_for method of XRTMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert xrt.is_datasource_for(xrt.data, xrt.meta)


def test_observatory(xrt):
    """Tests the observatory property of the XRTMap object."""
    assert xrt.observatory == "Hinode"


def test_measurement(xrt):
    """Tests the measurement property of the XRTMap object."""
    measurement = xrt.filter_wheel1_measurements[5].replace("_", " ")
    measurement += '-' + xrt.filter_wheel2_measurements[1].replace("_", " ")
    assert xrt.measurement == measurement


def test_wheel_measurements(xrt):
    """Tests the filter_wheel_measurements objects present
    in the XRTMap object."""
    assert (xrt.filter_wheel1_measurements ==
            ["Al_med", "Al_poly", "Be_med", "Be_thin", "C_poly", "Open"])
    assert (xrt.filter_wheel2_measurements ==
            ["Open", "Al_mesh", "Al_thick", "Be_thick", "Gband", "Ti_poly"])


def test_wcs(xrt):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        xrt.pixel_to_world(0*u.pix, 0*u.pix)
