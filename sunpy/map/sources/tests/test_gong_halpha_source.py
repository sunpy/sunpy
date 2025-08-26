import pytest
from numpy.testing import assert_equal

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources import GONGHalphaMap
from .helpers import _test_private_date_setters


@pytest.fixture
def gong_halpha():
    return get_dummy_map_from_header(get_test_filepath('gong_halpha.header'))


def test_fitstoGONGHAlphaMap(gong_halpha):
    """Tests the creation of GongSynopticMap using FITS."""
    assert isinstance(gong_halpha, GONGHalphaMap)


def test_is_datasource_for(gong_halpha):
    """Test the is_datasource_for method of GongSynopticMap."""
    assert gong_halpha.is_datasource_for(gong_halpha.data, gong_halpha.meta)


def test_observatory(gong_halpha):
    """Tests the observatory property of the GongSynopticMap object."""
    assert gong_halpha.observatory == "NSO-GONG"


def test_reference_date(gong_halpha):
    assert gong_halpha.reference_date.isot == "2024-02-16T00:00:02.000"


def test_date(gong_halpha):
    assert gong_halpha.date.isot == "2024-02-16T00:00:02.000"


def test_private_date_setters(gong_halpha):
    _test_private_date_setters(gong_halpha)


def test_scale(gong_halpha):
    """Tests the scale property of the GONGHalphaMap map."""
    assert_equal(gong_halpha.scale.axis1, 1.0794939681708389 * (u.arcsec / u.pix))
    assert_equal(gong_halpha.scale.axis2, 1.0794939681708389 * (u.arcsec / u.pix))


def test_rsun_obs(gong_halpha):
    assert gong_halpha.rsun_obs == 971.544571353755 * u.arcsec


def test_nickname(gong_halpha):
    """Tests the nickname property of the GONGHalphaMap map."""
    assert gong_halpha.nickname == "NSO-GONG, Big Bear, (USA)"


def test_earth_location(gong_halpha):
    assert u.isclose(gong_halpha._earth_location.lat, 34.26032998167749*u.deg, rtol=1e-15)
    assert u.isclose(gong_halpha._earth_location.lon, -116.92141999386104*u.deg, rtol=1e-15)


@pytest.mark.filterwarnings("ignore:Tried to get polar motions for times after IERS data is valid.")
@pytest.mark.filterwarnings("ignore:.*times are outside of range covered by IERS table.")
def test_observer_coordinate(gong_halpha):
    xyz_expected = [146712246479.363, -5563586.169750214, -17605285536.73928] * u.m
    assert u.isclose(gong_halpha.observer_coordinate.data.xyz, xyz_expected, rtol=1e-8).all()
