import pytest

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from sunpy.data.test import get_dummy_map_from_header_obj, get_header_from_file, get_test_filepath
from sunpy.map.sources.gong import GONGMagnetogramMap
from sunpy.time import parse_time
from .helpers import _test_private_date_setters


@pytest.fixture(params=["full-date-obs", "with-time-obs"])
def gong_magnetogram(request):
    header = get_header_from_file(get_test_filepath('gong_magnetogram.header'))
    if request.param == "full-date-obs":
        time_obs =  header.pop("TIME-OBS")
        date_obs = (parse_time(f"{header['DATE-OBS']} {time_obs}", scale='tai') + 19 *u.s).utc
        header["DATE-OBS"] = date_obs
    return get_dummy_map_from_header_obj(header)


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
    gong_magnetogram.wcs.pixel_to_world(0, 0)


def test_scale_units(gong_magnetogram):
    """Tests that the scale is reported in arcsec / pix, even though CDELTi are stored in spatial_units."""
    assert gong_magnetogram.scale.axis1.unit == u.arcsec / u.pix
    assert gong_magnetogram.scale.axis2.unit == u.arcsec / u.pix
    assert gong_magnetogram.spatial_units == (u.deg, u.deg)


@pytest.mark.parametrize("method", ["resample", "superpixel"])
def test_scale_after_resampling(gong_magnetogram, method):
    """Tests that halving the number of pixels doubles the scale and keeps the field of view."""
    ny, nx = gong_magnetogram.data.shape
    if method == "resample":
        smaller = gong_magnetogram.resample([nx // 2, ny // 2] * u.pix)
    else:
        smaller = gong_magnetogram.superpixel([2, 2] * u.pix)
    assert smaller.scale.axis1.unit == u.arcsec / u.pix
    assert smaller.scale.axis2.unit == u.arcsec / u.pix
    assert_quantity_allclose(smaller.scale.axis1, 2 * gong_magnetogram.scale.axis1)
    assert_quantity_allclose(smaller.scale.axis2, 2 * gong_magnetogram.scale.axis2)
    corner = gong_magnetogram.wcs.pixel_to_world(nx - 0.5, ny - 0.5)
    smaller_corner = smaller.wcs.pixel_to_world(nx / 2 - 0.5, ny / 2 - 0.5)
    assert_quantity_allclose(smaller_corner.Tx, corner.Tx)
    assert_quantity_allclose(smaller_corner.Ty, corner.Ty)


def test_scale_falls_back_to_limb_fit_without_cdelt(gong_magnetogram):
    """Tests that the scale is computed from the limb fit when CDELTi are removed."""
    limb_fit_scale = gong_magnetogram._limb_fit_scale
    del gong_magnetogram.meta['cdelt1']
    del gong_magnetogram.meta['cdelt2']
    assert_quantity_allclose(gong_magnetogram.scale.axis1, limb_fit_scale.axis1)
    assert_quantity_allclose(gong_magnetogram.scale.axis2, limb_fit_scale.axis2)


def test_scale_after_rotate_with_scale(gong_magnetogram):
    """Tests that rotating with a scale factor changes the scale, and keeps it in arcsec / pix."""
    rotated = gong_magnetogram.rotate(scale=2, missing=0)
    assert rotated.scale.axis1.unit == u.arcsec / u.pix
    assert rotated.scale.axis2.unit == u.arcsec / u.pix
    assert_quantity_allclose(rotated.scale.axis1, gong_magnetogram.scale.axis1 / 2)
    assert_quantity_allclose(rotated.scale.axis2, gong_magnetogram.scale.axis2 / 2)
