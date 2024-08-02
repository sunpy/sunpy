"""
Test cases for KCor Map subclass.
"""
import pytest

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.mlso import KCorMap
from .helpers import _test_private_date_setters


@pytest.fixture()
def kcor():
    """Creates an KCorMap from a FITS file."""
    return get_dummy_map_from_header(get_test_filepath("20181209_180305_kcor_l2.header"))


def test_kcormap_creation(kcor):
    """Tests the creation of KCorMap using FITS."""
    assert isinstance(kcor, KCorMap)


def test_is_datasource_for(kcor):
    """Test the is_datasource_for method of KCorMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert kcor.is_datasource_for(kcor.data, kcor.meta)


def test_measurement(kcor):
    """Tests the measurement property of the KCorMap object."""
    assert kcor.measurement == 735 * u.nm


def test_observatory(kcor):
    """Tests the observatory property of the KCorMap object."""
    assert kcor.observatory == "MLSO"


def test_norm_clip(kcor):
    # Tests that the default normalizer has clipping disabled
    assert not kcor.plot_settings['norm'].clip


def test_reference_date(kcor):
    assert kcor.reference_date.isot == "2018-12-09T18:03:05.000"


def test_date(kcor):
    assert kcor.date.isot == "2018-12-09T18:03:05.000"


def test_private_date_setters(kcor):
    _test_private_date_setters(kcor)


def test_wcs(kcor):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    kcor.pixel_to_world(0*u.pix, 0*u.pix)


def test_observer_coordinate(kcor):
    assert 'dsun_obs' not in kcor.meta

    # The test header triggers using the default observer coordinate
    observer = kcor.observer_coordinate.itrs.earth_location
    assert_quantity_allclose(observer.lon, kcor._earth_location.lon)
    assert_quantity_allclose(observer.lat, kcor._earth_location.lat)
    assert_quantity_allclose(observer.height, kcor._earth_location.height)

    # The test header will have a fully specified observer when DSUN_OBS is added
    kcor.meta['dsun_obs'] = (1*u.AU).to_value(u.m)
    kcor = KCorMap(kcor.data, kcor.meta)

    # The observer coordinate should now no longer be the default observer coordinate
    assert_quantity_allclose(kcor.observer_coordinate.radius, 1*u.AU)
    assert_quantity_allclose(kcor.observer_coordinate.lat, kcor.meta['crlt_obs']*u.deg)
