"""Tests for EUI Solar Orbiter Map"""

import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources import EUIMap
from .helpers import _test_private_date_setters

header_list = [
    get_test_filepath('solo_L1_eui-fsi304-image_20201021T145510206_V03.header'),
]


@pytest.fixture(scope="module", params=header_list)
def eui_map(request):
    return get_dummy_map_from_header(request.param)


def test_EUIMap(eui_map):
    assert isinstance(eui_map, EUIMap)


def test_reference_date(eui_map):
    assert eui_map.reference_date.isot == "2020-10-21T14:55:13.206"


def test_date(eui_map):
    assert eui_map.date.isot == "2020-10-21T14:55:10.206"


def test_private_date_setters(eui_map):
    _test_private_date_setters(eui_map)


def test_is_datasource_for(eui_map):
    assert eui_map.is_datasource_for(eui_map.data, eui_map.meta)


def test_observer_coordinate(eui_map):
    obs_coord = eui_map.observer_coordinate
    assert isinstance(obs_coord, SkyCoord)
    assert obs_coord.obstime.isot == eui_map.meta['date-avg']


def test_observatory(eui_map):
    assert eui_map.observatory == "Solar Orbiter"


def test_measurement(eui_map):
    assert eui_map.measurement == u.Quantity(304, 'angstrom')


def test_exposure_time(eui_map):
    assert eui_map.exposure_time == u.Quantity(6, 's')


def test_level_number(eui_map):
    assert eui_map.processing_level == 1


def test_unit(eui_map):
    assert eui_map.unit == u.Unit('DN')


def test_norm_clip(eui_map):
    # Tests that the default normalizer has clipping disabled
    assert not eui_map.plot_settings['norm'].clip


def test_wcs(eui_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    eui_map.pixel_to_world(0*u.pix, 0*u.pix)
