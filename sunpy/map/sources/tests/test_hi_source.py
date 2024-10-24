"""Test cases for STEREO Map subclasses.
This particular test file pertains to HIMap.
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.stereo import HIMap
from .helpers import _test_private_date_setters


@pytest.fixture
def hi_map():
    return get_dummy_map_from_header(get_test_filepath('hi_20110910_114721_s7h2A.header'))


def test_fitstoHI(hi_map):
    """Tests the creation of HIMap to fits"""
    assert isinstance(hi_map, HIMap)


def test_is_datasource_for(hi_map):
    """Test the is_data_source_for method of HIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert hi_map.is_datasource_for(hi_map.data, hi_map.meta)


def test_reference_date(hi_map):
    assert hi_map.reference_date.isot == "2011-09-10T11:47:46.004"


def test_date(hi_map):
    assert hi_map.date.isot == "2011-09-10T11:47:21.005"


def test_private_date_setters(hi_map):
    _test_private_date_setters(hi_map)


def test_measurement(hi_map):
    """Tests the measurement property of the HIMap object."""
    assert hi_map.measurement == "white-light"


def test_observatory(hi_map):
    """Tests the observatory property of the HIMap object."""
    assert hi_map.observatory == "STEREO A"


def test_norm_clip(hi_map):
    # Tests that the default normalizer has clipping disabled
    assert not hi_map.plot_settings['norm'].clip


def test_wcs(hi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    hi_map.pixel_to_world(0*u.pix, 0*u.pix)
