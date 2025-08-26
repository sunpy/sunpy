"""
Test cases for SUITMap subclass.
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.suit import SUITMap

__author__ = "Rahul Gopalakrishnan (rahulg.astro@gmail.com)"


@pytest.fixture
def suit_map():
    return get_dummy_map_from_header(get_test_filepath("SUT_T24_0847_000444_Lev1.0_2024-06-28T18.21.33.178_0971NB03.header"))


def test_suit_map(suit_map):
    """Tests the creation of SUITMap"""
    assert isinstance(suit_map, SUITMap)


def test_reference_date(suit_map):
    """ Tests the reference_date"""
    assert suit_map.reference_date.isot == "2024-06-28T18:21:33.178"


def test_date(suit_map):
    assert suit_map.date.isot == "2024-06-28T18:21:33.178"


def test_is_datasource_for(suit_map):
    """Tests the is_datasource_for method of SUITMap."""
    assert suit_map.is_datasource_for(suit_map.data, suit_map.meta)


def test_observatory(suit_map):
    """Tests the observatory property of the SUITMap object."""
    assert suit_map.observatory == "Aditya-L1"


def test_detector(suit_map):
    """Tests the detector property of the SUITMap object."""
    assert suit_map.detector == "SUIT"

def test_instrument(suit_map):
    """Tests the instrument property of the SUITMap object."""
    assert suit_map.instrument == "SUIT"

def test_norm_clip(suit_map):
    """Check for norm clipping"""
    assert not suit_map.plot_settings['norm'].clip


def test_wcs(suit_map):
    """Test WCS is valid and can transform from pixels to world coordinates"""
    suit_map.pixel_to_world(0*u.pix, 0*u.pix)
