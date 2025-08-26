"""
Test cases for SOHO EITL1Map subclass
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.mapbase import SpatialPair
from sunpy.map.sources.soho import EITL1Map
from .helpers import _test_private_date_setters


@pytest.fixture
def eit_l1_map():
    """Creates an EITL1Map"""
    return get_dummy_map_from_header(get_test_filepath("EIT_header/SOHO_EIT_171_20070601T120013_L1.header"))


def test_fitstoEIT(eit_l1_map):
    """Tests the creation of EITMap using FITS."""
    assert isinstance(eit_l1_map, EITL1Map)


def test_eitmap_coordinate_system(eit_l1_map):
    assert eit_l1_map.coordinate_system == SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')


def test_reference_date(eit_l1_map):
    assert eit_l1_map.reference_date.isot == "2007-06-01T11:59:05.180"


def test_date(eit_l1_map):
    assert eit_l1_map.date.isot == "2007-06-01T11:58:58.884"


def test_private_date_setters(eit_l1_map):
    _test_private_date_setters(eit_l1_map)


def test_is_datasource_for(eit_l1_map):
    """Test the is_datasource_for method of EITMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert eit_l1_map.is_datasource_for(eit_l1_map.data, eit_l1_map.meta)


def test_observatory(eit_l1_map):
    """Tests the observatory property of the EITMap object."""
    assert eit_l1_map.observatory == "Solar and Heliospheric Observatory (SOHO)"

def test_instrument(eit_l1_map):
    """Tests the instrument property of the EITMap object."""
    assert eit_l1_map.instrument == "Extreme-ultraviolet Imaging Telescope (EIT)"

def test_measurement(eit_l1_map):
    """Tests the measurement property of the EITMap object."""
    assert eit_l1_map.measurement == 171*u.angstrom


def test_rsun(eit_l1_map):
    """Tests the measurement property of the EITMap object."""
    assert u.allclose(eit_l1_map.rsun_obs, 953.98737017*u.arcsec)


def test_norm_clip(eit_l1_map):
    # Tests that the default normalizer has clipping disabled
    assert not eit_l1_map.plot_settings['norm'].clip


def test_wcs(eit_l1_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    eit_l1_map.pixel_to_world(0*u.pix, 0*u.pix)


def test_l1_eit_l1_map_colormaps():
    eit_l1_map = get_dummy_map_from_header(get_test_filepath("EIT_header/SOHO_EIT_171_20070601T120013_L1.header"))
    assert eit_l1_map.plot_settings["cmap"] == "sohoeit171"
    eit_l1_map = get_dummy_map_from_header(get_test_filepath("EIT_header/SOHO_EIT_195_20070601T121346_L1.header"))
    assert eit_l1_map.plot_settings["cmap"] == "sohoeit195"
    eit_l1_map = get_dummy_map_from_header(get_test_filepath("EIT_header/SOHO_EIT_284_20070601T120607_L1.header"))
    assert eit_l1_map.plot_settings["cmap"] == "sohoeit284"
    eit_l1_map = get_dummy_map_from_header(get_test_filepath("EIT_header/SOHO_EIT_304_20070601T121937_L1.header"))
    assert eit_l1_map.plot_settings["cmap"] == "sohoeit304"
