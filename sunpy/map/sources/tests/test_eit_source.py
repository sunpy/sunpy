"""
Test cases for SOHO EITMap subclass
"""
import pytest

import astropy.units as u

import sunpy.map
from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.mapbase import SpatialPair
from sunpy.map.sources.soho import EITMap
from .helpers import _test_private_date_setters

__author__ = "Pritish C. (VaticanCameos)"


@pytest.fixture()
def eit_map():
    """Creates an EITMap from a FITS file."""
    return sunpy.map.Map(get_test_filepath("EIT/efz20040301.000010_s.fits"))


def test_fitstoEIT(eit_map):
    """Tests the creation of EITMap using FITS."""
    assert isinstance(eit_map, EITMap)


def test_eitmap_coordinate_system(eit_map):
    assert eit_map.coordinate_system ==  SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')


def test_reference_date(eit_map):
    assert eit_map.reference_date.isot == "2004-03-01T00:00:10.515"


def test_date(eit_map):
    assert eit_map.date.isot == "2004-03-01T00:00:10.515"


def test_private_date_setters(eit_map):
    _test_private_date_setters(eit_map)


def test_is_datasource_for(eit_map):
    """Test the is_datasource_for method of EITMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert eit_map.is_datasource_for(eit_map.data, eit_map.meta)


def test_observatory(eit_map):
    """Tests the observatory property of the EITMap object."""
    assert eit_map.observatory == "SOHO"


def test_measurement(eit_map):
    """Tests the measurement property of the EITMap object."""
    assert eit_map.measurement.value in [195, 171]


def test_rsun(eit_map):
    """Tests the measurement property of the EITMap object."""
    assert u.allclose(eit_map.rsun_obs, 979.0701*u.arcsec)


def test_norm_clip(eit_map):
    # Tests that the default normalizer has clipping disabled
    assert not eit_map.plot_settings['norm'].clip


def test_wcs(eit_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    eit_map.pixel_to_world(0*u.pix, 0*u.pix)


def test_old_eit_date():
    eit_map = get_dummy_map_from_header(get_test_filepath("seit_00171_fd_19961211_1900.header"))
    assert eit_map.date.value == '1996-12-11T19:00:14.254'
    assert eit_map.reference_date.value == '1996-12-11T19:00:14.254'
