"""Tests for PHI Solar Orbiter Map"""

import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources import PHIMap
from .helpers import _test_private_date_setters

header_list = [
    get_test_filepath('solo_L2_phi-hrt-blos_20241004T003104_V202506050052_0450040601.header'),
    get_test_filepath('solo_LL02_phi-fdt-blos_20240305T041509_V202405151730C_0403057611.header')
]


@pytest.fixture(scope="module", params=header_list)
def phi_map(request):
    return get_dummy_map_from_header(request.param)


def test_PHIMap(phi_map):
    assert isinstance(phi_map, PHIMap)


def test_reference_date(phi_map):
    assert phi_map.reference_date.isot == "2020-10-21T14:55:13.206"


def test_date(phi_map):
    assert phi_map.date.isot == "2020-10-21T14:55:10.206"


def test_private_date_setters(phi_map):
    _test_private_date_setters(phi_map)


def test_is_datasource_for(phi_map):
    assert phi_map.is_datasource_for(phi_map.data, phi_map.meta)


def test_observer_coordinate(phi_map):
    obs_coord = phi_map.observer_coordinate
    assert isinstance(obs_coord, SkyCoord)
    assert obs_coord.obstime.isot == phi_map.meta['date-avg']


def test_observatory(phi_map):
    assert phi_map.observatory == "Solar Orbiter"


def test_measurement(phi_map):
    assert phi_map.measurement == u.Quantity(304, 'angstrom')


def test_exposure_time(phi_map):
    assert phi_map.exposure_time == u.Quantity(6, 's')


def test_level_number(phi_map):
    assert phi_map.processing_level == 1


def test_unit(phi_map):
    assert phi_map.unit == u.Unit('DN')


def test_norm_clip(phi_map):
    # Tests that the default normalizer has clipping disabled
    assert not phi_map.plot_settings['norm'].clip


def test_wcs(phi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    phi_map.pixel_to_world(0*u.pix, 0*u.pix)
