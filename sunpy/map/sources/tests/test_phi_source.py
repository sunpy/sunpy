"""Tests for PHI Solar Orbiter Map"""

import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
import warnings

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources import PHIMap
from sunpy.util.exceptions import SunpyUserWarning
from .helpers import _test_private_date_setters

hrt_header_list = [
    get_test_filepath('solo_L2_phi-hrt-vlos_20241004T003104_V202506050052_0450040601.header'),
    get_test_filepath('solo_L2_phi-hrt-blos_20241004T003104_V202506050052_0450040601.header'),
    get_test_filepath('solo_L2_phi-hrt-bmag_20241004T003104_V202506050052_0450040601.header'),
    get_test_filepath('solo_L2_phi-hrt-binc_20241004T003104_V202506050052_0450040601.header'),
    get_test_filepath('solo_L2_phi-hrt-bazi_20241004T003104_V202506050052_0450040601.header'),
    get_test_filepath('solo_L2_phi-hrt-icnt_20241004T003104_V202506050052_0450040601.header'),
    get_test_filepath('solo_L2_phi-hrt-vlos_20220307T000009_V202208311927_0243070101.header'),
    get_test_filepath('solo_L2_phi-hrt-blos_20220307T000009_V202208311927_0243070101.header'),
    get_test_filepath('solo_L2_phi-hrt-bmag_20220307T000009_V202208311927_0243070101.header'),
    get_test_filepath('solo_L2_phi-hrt-binc_20220307T000009_V202208311927_0243070101.header'),
    get_test_filepath('solo_L2_phi-hrt-bazi_20220307T000009_V202208311927_0243070101.header'),
    get_test_filepath('solo_L2_phi-hrt-icnt_20220307T000009_V202208311927_0243070101.header'),
    
]

expected_hrt_refdates_list = 6*['2024-10-04T00:31:45.499'] + 6*['2022-03-07T00:00:32.393']

expected_hrt_dates_list = 6*['2024-10-04T00:31:04.322'] + 6*['2022-03-07T00:00:09.388']

expected_hrt_measurement_list = [
    'VLOS',
    'BLOS',
    'BMAG',
    'BINC',
    'BAZI',
    'ICNT',
    'VLOS',
    'BLOS',
    'BMAG',
    'BINC',
    'BAZI',
    'ICNT',
]

expected_hrt_units_list = [
    u.Unit('km/s'),
    u.Unit('G'),
    u.Unit('G'),
    u.Unit('deg'),
    u.Unit('deg'),
    u.Unit(''),
    u.Unit('km/s'),
    u.Unit('G'),
    u.Unit('G'),
    u.Unit('deg'),
    u.Unit('deg'),
    u.Unit(''),
]

test_hrt_cal_wcs_warning_header_list = [
    get_test_filepath('solo_L2_phi-hrt-vlos_20220307T000009_V202208311927_0243070101.header')
]

stokes_header_list = [
    get_test_filepath('solo_L2_phi-hrt-stokes_20241004T003104_V202506050052_0450040601.header'),
    get_test_filepath('solo_L2_phi-hrt-stokes_20220307T000009_V202208311927_0243070101.header'),
]

fdt_header_list = [
    get_test_filepath('solo_LL02_phi-fdt-blos_20240305T041509_V202405151730C_0403057611.header'),
]


@pytest.fixture(scope="module", params=hrt_header_list)
def phi_map_hrt(request):
    import warnings
    from sunpy.util.exceptions import SunpyUserWarning
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=SunpyUserWarning,
            message=r".*WCS of this SO/PHI-HRT PHIMap.*calibrated.*"
        )
        warnings.filterwarnings(
            "ignore",
            category=SunpyUserWarning,
            message=r".*may not be fully calibrated.*"
        )
        return get_dummy_map_from_header(request.param)


@pytest.fixture(scope="module", params=stokes_header_list)
def phi_map_stokes(request):
    return get_dummy_map_from_header(request.param)


@pytest.fixture(scope="module", params=test_hrt_cal_wcs_warning_header_list)
def test_phi_map_hrt_wcs_warning(request):
    with pytest.warns(SunpyUserWarning):
        get_dummy_map_from_header(request.param)


def test_hrt_PHIMap(phi_map_hrt):
    assert isinstance(phi_map_hrt, PHIMap)


def test_stokes_PHIMap(phi_map_stokes):
    assert not isinstance(phi_map_stokes, PHIMap)


@pytest.mark.parametrize(
        ('phi_map_hrt', 'expected_refdate'),
        list(zip(hrt_header_list, expected_hrt_refdates_list)),
        indirect=['phi_map_hrt']
        )
def test_reference_date(phi_map_hrt, expected_refdate):
    assert phi_map_hrt.reference_date.isot == expected_refdate


@pytest.mark.parametrize(
        ('phi_map_hrt', 'expected_date'),
        list(zip(hrt_header_list, expected_hrt_dates_list)),
        indirect=['phi_map_hrt']
        )
def test_date(phi_map_hrt, expected_date):
    assert phi_map_hrt.date.isot == expected_date


def test_private_date_setters(phi_hrt_map):
    _test_private_date_setters(phi_hrt_map)


@pytest.mark.parametrize(
        ('phi_map_hrt', 'expected_measurement'),
        list(zip(hrt_header_list, expected_hrt_measurement_list)),
        indirect=['phi_map_hrt']
        )
def test_measurement(phi_map_hrt, expected_measurement):
    assert phi_map_hrt.measurement == expected_measurement


def test_is_datasource_for(phi_map_hrt):
    assert phi_map_hrt.is_datasource_for(phi_map_hrt.data, phi_map_hrt.meta)


def test_observer_coordinate(phi_map_hrt):
    obs_coord = phi_map_hrt.observer_coordinate
    assert isinstance(obs_coord, SkyCoord)
    assert obs_coord.obstime.isot == phi_map_hrt.meta['date-avg']


def test_observatory(phi_map_hrt):
    assert phi_map_hrt.observatory == "Solar Orbiter"


def test_level_number(phi_map_hrt):
    assert phi_map_hrt.processing_level == 2


@pytest.mark.parametrize(
        ('phi_map_hrt', 'expected_unit'),
        list(zip(hrt_header_list, expected_hrt_units_list)),
        indirect=['phi_map_hrt']
        )
def test_unit(phi_map_hrt, expected_unit):
    assert phi_map_hrt.unit == expected_unit


def test_wcs(phi_map_hrt):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    phi_map_hrt.pixel_to_world(0*u.pix, 0*u.pix)
