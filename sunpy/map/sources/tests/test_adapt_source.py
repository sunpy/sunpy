import pytest

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.gong import ADAPTMap


@pytest.fixture
def adapt_map():
    return get_dummy_map_from_header(get_test_filepath('adapt.header'))


def test_fitstoADAPTMap(adapt_map):
    """Tests the creation of ADAPTMap using FITS."""
    assert isinstance(adapt_map, ADAPTMap)


def test_is_datasource_for(adapt_map):
    """Test the is_datasource_for method of GongSynopticMap."""
    assert adapt_map.is_datasource_for(adapt_map.data, adapt_map.meta)


def test_date(adapt_map):
    """Check that accessing the date doesn't raise a warning."""
    adapt_map.date


def test_date_maptime(adapt_map):
    assert adapt_map.meta['date-obs'] == adapt_map.meta['maptime']


def test_ctype(adapt_map):
    assert adapt_map.meta['ctype1'] == 'CRLN-CAR'
    assert adapt_map.meta['ctype2'] == 'CRLT-CAR'
