import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.hinode import SOTMap, XRTMap
from sunpy.util.exceptions import SunpyMetadataWarning


@pytest.fixture
def xrt_map():
    return get_dummy_map_from_header(get_test_filepath("HinodeXRT.header"))


@pytest.fixture
def sot_map():
    return get_dummy_map_from_header(get_test_filepath("HinodeSOT.header"))


def test_xrt_map(xrt_map):
    assert isinstance(xrt_map, XRTMap)


def test_is_datasource_for_xrt(xrt_map):
    assert xrt_map.is_datasource_for(xrt_map.data, xrt_map.meta)


def test_observatory_xrt(xrt_map):
    assert xrt_map.observatory == "Hinode"


def test_measurement_xrt(xrt_map):
    assert xrt_map.measurement == "Be thin-Open"


def test_wheel_measurements_xrt(xrt_map):
    assert (xrt_map.filter_wheel1_measurements ==
            ["Al_med", "Al_poly", "Be_med", "Be_thin", "C_poly", "Open"])
    assert (xrt_map.filter_wheel2_measurements ==
            ["Open", "Al_mesh", "Al_thick", "Be_thick", "Gband", "Ti_poly"])


def test_wcs_xrt(xrt_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    xrt_map.pixel_to_world(0*u.pix, 0*u.pix)


def test_sot_map(sot_map):
    assert isinstance(sot_map, SOTMap)


def test_is_datasource_for_sot(sot_map):
    assert sot_map.is_datasource_for(sot_map.data, sot_map.meta)


def test_observatory_sot(sot_map):
    assert sot_map.observatory == "Hinode"


def test_measurement_sot(sot_map):
    assert sot_map.measurement is None


def test_instruments_sot(sot_map):
    assert sot_map.instruments == ['SOT/WB',
                                   'SOT/NB', 'SOT/SP', 'SOT/CT']


def test_waves_sot(sot_map):
    assert sot_map.waves == ['6302A', 'BFI no move',
                             'CN bandhead 3883', 'Ca II H line',
                             'G band 4305', 'NFI no move', 'TF Fe I 6302',
                             'TF Mg I 5172', 'TF Na I 5896',
                             'blue cont 4504', 'green cont 5550',
                             'red cont 6684']


def test_obstype_sot(sot_map):
    assert sot_map.observation_type == ['FG (simple)',
                                        'FG focus scan', 'FG shuttered I and V',
                                        'FG shutterless I and V', 'FG shutterless I and V with 0.2s intervals',
                                        'FG shutterless Stokes', 'SP IQUV 4D array']


def test_wcs_sot(sot_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        sot_map.pixel_to_world(0*u.pix, 0*u.pix)
