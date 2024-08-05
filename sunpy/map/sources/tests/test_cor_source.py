"""
Test cases for STEREO CORMap subclass.
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.stereo import CORMap
from .helpers import _test_private_date_setters

__author__ = 'Pritish C. (VaticanCameos)'


@pytest.fixture
def cor_map():
    header_file = get_test_filepath("cor1_20090615_000500_s4c1A.header")
    return get_dummy_map_from_header(header_file)


def test_fits_to_cor(cor_map):
    """Tests the creation of CORMap using FITS."""
    assert isinstance(cor_map, CORMap)


def test_reference_date(cor_map):
    assert cor_map.reference_date.isot == "2009-06-15T00:05:00.855"


def test_date(cor_map):
    assert cor_map.reference_date.isot == "2009-06-15T00:05:00.855"


def test_private_date_setters(cor_map):
    _test_private_date_setters(cor_map)


def test_is_datasource_for(cor_map):
    """Test the is_datasource_for method of CORMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert cor_map.is_datasource_for(cor_map.data, cor_map.meta)


def test_measurement(cor_map):
    """Tests the measurement property of the CORMap object."""
    assert cor_map.measurement == "white-light"


def test_observatory(cor_map):
    """Tests the observatory property of the CORMap object."""
    assert cor_map.observatory == "STEREO A"


def test_norm_clip(cor_map):
    # Tests that the default normalizer has clipping disabled
    assert not cor_map.plot_settings['norm'].clip


def test_wcs(cor_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    cor_map.pixel_to_world(0*u.pix, 0*u.pix)
