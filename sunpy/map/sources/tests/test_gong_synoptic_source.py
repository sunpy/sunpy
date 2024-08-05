import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.gong import GONGSynopticMap
from .helpers import _test_private_date_setters


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
    assert gong_synoptic.measurement == 676.8 * u.nm


def test_reference_date(gong_synoptic):
    assert gong_synoptic.reference_date.isot == "2023-09-30T06:44:00.000"


def test_date(gong_synoptic):
    assert gong_synoptic.date.isot == "2023-09-30T06:44:00.000"


def test_private_date_setters(gong_synoptic):
    _test_private_date_setters(gong_synoptic)


def test_unit(gong_synoptic):
    assert gong_synoptic.unit == u.G
    assert gong_synoptic.unit == u.Unit("Mx/cm^2")
    assert gong_synoptic.unit.to_string() == 'G'


def test_spatial_units(gong_synoptic):
    assert gong_synoptic.spatial_units[0] == u.deg
    assert gong_synoptic.spatial_units[1] == u.deg


def test_wcs(gong_synoptic):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    gong_synoptic.pixel_to_world(0*u.pix, 0*u.pix)
