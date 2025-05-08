import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.gong import GONGMagnetogramMap
from .helpers import _test_private_date_setters


@pytest.fixture
def gong_magnetogram():
    return get_dummy_map_from_header(get_test_filepath('gong_magnetogram.header'))


def test_fitstoGONGSynoptic(gong_magnetogram):
    """Tests the creation of GONGMagnetogramMap using FITS."""
    assert isinstance(gong_magnetogram, GONGMagnetogramMap)


def test_is_datasource_for(gong_magnetogram):
    """Test the is_datasource_for method of GONGMagnetogramMap."""
    assert gong_magnetogram.is_datasource_for(gong_magnetogram.data, gong_magnetogram.meta)


def test_observatory(gong_magnetogram):
    """Tests the observatory property of the GONGMagnetogramMap object."""
    assert gong_magnetogram.observatory == "NSO-GONG"


def test_measurement(gong_magnetogram):
    """Tests the measurement property of the GONGMagnetogramMap object."""
    assert gong_magnetogram.measurement == 676.8 * u.nm


def test_reference_date(gong_magnetogram):
    # 2025-01-09T01:51:57
    assert gong_magnetogram.reference_date.isot == "2025-01-09T01:34:58.008"


def test_date(gong_magnetogram):
    # The time in the GONG header is in a ISO format but in GPS scale
    # GPS time (date-obs) 2025-01-09T01:51:57 should be UTC 2025-01-09T01:34:58.008 (-18 seconds)
    assert gong_magnetogram.date.isot == "2025-01-09T01:34:58.008"


def test_nickname(gong_magnetogram):
    assert gong_magnetogram.nickname == "NSO-GONG, Learmonth, (AUS)"


def test_private_date_setters(gong_magnetogram):
    _test_private_date_setters(gong_magnetogram)


def test_unit(gong_magnetogram):
    assert gong_magnetogram.unit == u.G
    assert gong_magnetogram.unit == u.Unit("Mx/cm^2")
    assert gong_magnetogram.unit.to_string() == 'G'


def test_spatial_units(gong_magnetogram):
    assert gong_magnetogram.spatial_units[0] == u.deg
    assert gong_magnetogram.spatial_units[1] == u.deg


def test_rsun_obs(gong_magnetogram):
    assert u.allclose(gong_magnetogram.rsun_obs, 975.86818683 * u.arcsec)


def test_wcs(gong_magnetogram):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    gong_magnetogram.pixel_to_world(0*u.pix, 0*u.pix)
