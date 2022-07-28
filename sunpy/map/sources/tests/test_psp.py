import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources import WISPRMap


@pytest.fixture
def wispr_map(request):
    return get_dummy_map_from_header(get_test_filepath("solo_eui_fsi_304_wispr.header"))


def test_wispr_map(wispr_map):
    assert isinstance(wispr_map, WISPRMap)


def test_is_datasource_for_wispr(wispr_map):
    assert wispr_map.is_datasource_for(wispr_map.data, wispr_map.meta)


def test_observer_coordinate_wispr(wispr_map):
    obs_coord = wispr_map.observer_coordinate
    assert isinstance(obs_coord, SkyCoord)
    assert obs_coord.obstime.isot == wispr_map.meta['date-obs']


def test_observatory_wispr(wispr_map):
    assert wispr_map.observatory == "Parker Solar Probe"


def test_measurement_wispr(wispr_map):
    assert wispr_map.measurement is None


def test_wavelength_wispr(wispr_map):
    assert wispr_map.wavelength is None


def test_exposure_time_wispr(wispr_map):
    assert wispr_map.exposure_time == u.Quantity(700, 's')


def test_level_number_wispr(wispr_map):
    assert wispr_map.processing_level == 1


def test_detector_wispr(wispr_map):
    assert wispr_map.detector == 2


def test_unit_wispr(wispr_map):
    assert wispr_map.unit == u.Unit('ct')


def test_norm_clip_wispr(wispr_map):
    # Tests that the default normalizer has clipping disabled
    assert not wispr_map.plot_settings['norm'].clip


def test_name_wispr(wispr_map):
    assert wispr_map.name == 'WISPR 2 2020-01-25 00:02:29'


def test_wcs_wispr(wispr_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    wispr_map.pixel_to_world(0*u.pix, 0*u.pix)
