import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.gong import GONGSynopticMap


@pytest.fixture
def gong_synoptic():
    return get_dummy_map_from_header(get_test_filepath('gong_synoptic.header'))


def test_fitstoGONGSynoptic(gong_synoptic):
    """Tests the creation of GongSynopticMap using FITS."""
    assert isinstance(gong_synoptic, GONGSynopticMap)


def test_is_datasource_for(gong_synoptic):
    """Test the is_datasource_for method of GongSynopticMap."""
    assert gong_synoptic.is_datasource_for(gong_synoptic.data, gong_synoptic.meta)


def test_observatory(gong_synoptic):
    """Tests the observatory property of the GongSynopticMap object."""
    assert gong_synoptic.observatory == "NSO-GONG"


def test_measurement(gong_synoptic):
    """Tests the measurement property of the GongSynopticMap object."""
    assert gong_synoptic.measurement == 676.8


def test_date(gong_synoptic):
    """Check that accessing the date doesn't raise a warning."""
    gong_synoptic.date


def test_unit(gong_synoptic):
    assert gong_synoptic.unit == u.G
    assert gong_synoptic.unit == u.Unit("Mx/cm^2")
    assert gong_synoptic.unit.to_string() == 'G'
    gong_synoptic.meta['bunit'] = 'm'
    assert gong_synoptic.unit == u.m


def test_coordinate_system(gong_synoptic):
    assert gong_synoptic.meta['CUNIT1'] == 'deg'
    assert gong_synoptic.meta['CUNIT2'] == 'deg'


def test_observer(gong_synoptic):
    assert gong_synoptic.meta['hglt_obs'] is not None
    assert gong_synoptic.meta['hgln_obs'] is not None
    assert gong_synoptic.meta['dsun_obs'] is not None
