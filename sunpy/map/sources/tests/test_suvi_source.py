"""
Test cases for SUVI Map subclass.
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.suvi import SUVIMap
from .helpers import _test_private_date_setters


@pytest.fixture()
def suvi():
    """Creates an SUVIMap from a FITS file."""
    path = get_test_filepath("dr_suvi-l2-ci195_g16_s20190403T093200Z_e20190403T093600Z_v1-0-0_rebinned.header")
    return get_dummy_map_from_header(path)


def test_suvimap_creation(suvi):
    """Tests the creation of SUVIMap using FITS."""
    assert isinstance(suvi, SUVIMap)


def test_reference_date(suvi):
    assert suvi.reference_date.isot == "2019-04-03T09:32:33.340"


def test_date(suvi):
    assert suvi.date.isot == "2019-04-03T09:32:33.340"


def test_private_date_setters(suvi):
    _test_private_date_setters(suvi)


def test_is_datasource_for(suvi):
    """Test the is_datasource_for method of SUVIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert suvi.is_datasource_for(suvi.data, suvi.meta)


def test_observatory(suvi):
    """Tests the observatory property of the SUVIMap object."""
    assert suvi.observatory == "GOES-R"


def test_detector(suvi):
    """Tests the detector property of the SUVIMap object."""
    assert suvi.detector == "SUVI"


def test_norm_clip(suvi):
    # Tests that the default normalizer has clipping disabled
    assert not suvi.plot_settings['norm'].clip


# SUVI provides observer coordinate information in an OBSGEO system, so this test
# needs remote data to access the latest IERS table to do a coordinate transformation from
# OBSGEO to heliographic Stonyhurst coordinates.
@pytest.mark.remote_data
def test_wcs(suvi):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    suvi.pixel_to_world(0*u.pix, 0*u.pix)
