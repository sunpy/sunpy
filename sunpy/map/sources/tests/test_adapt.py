
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.adapt import ADAPTMap
from sunpy.util.exceptions import SunpyMetadataWarning


@pytest.fixture
def adapt_map():
    header_file = get_test_filepath("adapt.header")
    return get_dummy_map_from_header(header_file)


def test_fitstoadapt(adapt_map):
    """Tests the creation of ADAPTMap using FITS."""
    assert isinstance(adapt_map, ADAPTMap)


def test_is_datasource_for(adapt_map):
    """Test the is_datasource_for method of ADAPTMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert adapt_map.is_datasource_for(adapt_map.data, adapt_map.meta)


def test_measurement(adapt_map):
    """Tests the measurement property of the ADAPTMap object."""
    assert adapt_map.measurement is None


def test_observatory(adapt_map):
    """Tests the observatory property of the ADAPTMap object."""
    assert adapt_map.observatory == ""


def test_wcs(adapt_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='assuming Earth-based observer'):
        adapt_map.pixel_to_world(0*u.pix, 0*u.pix)
