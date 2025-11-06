"""
Test cases for SDO HMIMap subclass.
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map import Map
from sunpy.map.sources.sdo import HMIMap
from .helpers import _test_private_date_setters

__author__ = 'Pritish C. (VaticanCameos)'


@pytest.fixture
def hmi_map():
    return Map(get_test_filepath('resampled_hmi.fits'))


@pytest.fixture
def hmi_bharp_map():
    return get_dummy_map_from_header(get_test_filepath('hmi_bharp_vlos_mag.header'))


@pytest.fixture
def hmi_cea_sharp_map():
    return get_dummy_map_from_header(get_test_filepath('hmi_cea_sharp_magnetogram.header'))


@pytest.fixture
def hmi_sharp_map():
    return get_dummy_map_from_header(get_test_filepath('hmi_sharp_magnetogram.header'))


def test_fitstoHMI(hmi_map, hmi_bharp_map, hmi_cea_sharp_map, hmi_sharp_map):
    """Tests the creation of HMIMap using FITS."""
    assert isinstance(hmi_map, HMIMap)
    assert isinstance(hmi_bharp_map, HMIMap)
    assert isinstance(hmi_cea_sharp_map, HMIMap)
    assert isinstance(hmi_sharp_map, HMIMap)


def test_is_datasource_for(hmi_map, hmi_bharp_map, hmi_cea_sharp_map, hmi_sharp_map):
    """Test the is_datasource_for method of HMIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object, which in this case is
    hmi.meta."""
    assert hmi_map.is_datasource_for(hmi_map.data, hmi_map.meta)
    assert hmi_bharp_map.is_datasource_for(hmi_bharp_map.data, hmi_bharp_map.meta)
    assert hmi_cea_sharp_map.is_datasource_for(hmi_cea_sharp_map.data, hmi_cea_sharp_map.meta)
    assert hmi_sharp_map.is_datasource_for(hmi_sharp_map.data, hmi_sharp_map.meta)


def test_reference_date(hmi_map, hmi_bharp_map, hmi_cea_sharp_map, hmi_sharp_map):
    assert hmi_map.reference_date.isot == "2014-03-01T00:01:25.000"
    assert hmi_bharp_map.reference_date.isot == "2014-06-09T23:48:07.532"
    assert hmi_cea_sharp_map.reference_date.isot == "2024-06-28T00:00:08.212"
    assert hmi_sharp_map.reference_date.isot == "2024-06-28T00:00:08.212"


def test_date(hmi_map, hmi_bharp_map, hmi_cea_sharp_map, hmi_sharp_map):
    assert hmi_map.date.isot == "2014-03-01T00:00:27.900"
    assert hmi_bharp_map.date.isot == "2014-06-09T23:46:25.000"
    assert hmi_cea_sharp_map.date.isot == "2024-06-27T23:58:46.200"
    assert hmi_sharp_map.date.isot == "2024-06-27T23:58:46.200"


def test_private_date_setters(hmi_map, hmi_bharp_map, hmi_cea_sharp_map, hmi_sharp_map):
    _test_private_date_setters(hmi_map)
    _test_private_date_setters(hmi_bharp_map)
    _test_private_date_setters(hmi_cea_sharp_map)
    _test_private_date_setters(hmi_sharp_map)


def test_observatory(hmi_map, hmi_bharp_map, hmi_cea_sharp_map, hmi_sharp_map):
    """Tests the observatory property of the HMIMap object."""
    assert hmi_map.observatory == "SDO"
    assert hmi_bharp_map.observatory == "SDO"
    assert hmi_cea_sharp_map.observatory == "SDO"
    assert hmi_sharp_map.observatory == "SDO"


def test_measurement(hmi_map, hmi_bharp_map, hmi_cea_sharp_map, hmi_sharp_map):
    """Tests the measurement property of the HMIMap object."""
    assert hmi_map.measurement == "continuum"
    assert hmi_bharp_map.measurement == "observable"
    assert hmi_cea_sharp_map.measurement == "observable"
    assert hmi_sharp_map.measurement == "observable"


def test_wavelength(hmi_map, hmi_bharp_map, hmi_cea_sharp_map, hmi_sharp_map):
    assert hmi_map.wavelength == 6173  # The test file actually has an "" WAVEUNIT key
    assert hmi_bharp_map.wavelength == 6173 * u.AA
    assert hmi_cea_sharp_map.wavelength == 6173 * u.AA
    assert hmi_sharp_map.wavelength == 6173 * u.AA


def test_unit(hmi_cea_sharp_map, hmi_sharp_map):
    """Test that HMI source handles units not in FITS standard"""
    assert hmi_cea_sharp_map.unit == u.Unit('Mx/cm2')
    assert hmi_sharp_map.unit == u.Unit('Mx/cm2')


def test_wcs(hmi_map, hmi_bharp_map, hmi_cea_sharp_map, hmi_sharp_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    hmi_map.pixel_to_world(0*u.pix, 0*u.pix)
    hmi_bharp_map.pixel_to_world(0*u.pix, 0*u.pix)
    hmi_cea_sharp_map.pixel_to_world(0*u.pix, 0*u.pix)
    hmi_sharp_map.pixel_to_world(0*u.pix, 0*u.pix)
