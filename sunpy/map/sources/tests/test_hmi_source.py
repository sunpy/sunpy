"""
Test cases for SDO HMIMap subclass.
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_test_filepath
from sunpy.map import Map
from sunpy.map.sources.sdo import HMIMap

__author__ = 'Pritish C. (VaticanCameos)'


@pytest.fixture
def hmi_map():
    return Map(get_test_filepath('resampled_hmi.fits'))


def test_fitstoHMI(hmi_map):
    """Tests the creation of HMIMap using FITS."""
    assert isinstance(hmi_map, HMIMap)


def test_is_datasource_for(hmi_map):
    """Test the is_datasource_for method of HMIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object, which in this case is
    hmi.meta."""
    assert hmi_map.is_datasource_for(hmi_map.data, hmi_map.meta)


def test_observatory(hmi_map):
    """Tests the observatory property of the HMIMap object."""
    assert hmi_map.observatory == "SDO"


def test_measurement(hmi_map):
    """Tests the measurement property of the HMIMap object."""
    assert hmi_map.measurement == "continuum"


def test_wcs(hmi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    hmi_map.pixel_to_world(0*u.pix, 0*u.pix)
