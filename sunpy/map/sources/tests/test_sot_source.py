
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.mapbase import SpatialPair
from sunpy.map.sources.hinode import SOTMap
from sunpy.util.exceptions import SunpyMetadataWarning
from .helpers import _test_private_date_setters


@pytest.fixture
def sot_fg_wb_map():
    return get_dummy_map_from_header(get_test_filepath("HinodeSOT/HinodeSOT_FG_WB_G_band_4305_L0.header"))


SOT_LOADABLE_HEADERS = [
    ("HinodeSOT/HinodeSOT_FG_WB_G_band_4305_L0.header",
     "SOT/WB", "G band 4305"),
    ("HinodeSOT/HinodeSOT_FG_WB_CaIIH_L1.header",
     "SOT/WB", "Ca II H line"),
    ("HinodeSOT/HinodeSOT_FG_NB_TF_Na_I_5896_L1.header",
     "SOT/NB", "TF Na I 5896"),
]


@pytest.fixture(params=SOT_LOADABLE_HEADERS,
                ids=[h[0] for h in SOT_LOADABLE_HEADERS])
def sot_map(request):
    header_file, _, _ = request.param
    m = get_dummy_map_from_header(get_test_filepath(header_file))
    return m, request.param


def test_sot_map_is_sotmap(sot_map):
    m, _ = sot_map
    assert isinstance(m, SOTMap)


def test_sot_map_observatory(sot_map):
    m, _ = sot_map
    assert m.observatory == "Hinode"


def test_sot_map_detector(sot_map):
    m, _ = sot_map
    assert m.detector == "SOT"


def test_sot_map_instrument(sot_map):
    m, (_, expected_instrume, _) = sot_map
    assert m.instrument == expected_instrume


def test_sot_map_wave_in_waves(sot_map):
    m, (_, _, expected_wave) = sot_map
    assert expected_wave in SOTMap.Waves


def test_sot_map_coordinate_system(sot_map):
    m, _ = sot_map
    assert m.coordinate_system == SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')


def test_fitstoSOT(sot_fg_wb_map):
    assert isinstance(sot_fg_wb_map, SOTMap)


def test_sot_coordinate_system(sot_fg_wb_map):
    assert sot_fg_wb_map.coordinate_system == SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')


def test_is_datasource_for(sot_fg_wb_map):
    assert sot_fg_wb_map.is_datasource_for(sot_fg_wb_map.data, sot_fg_wb_map.meta)


def test_observatory(sot_fg_wb_map):
    assert sot_fg_wb_map.observatory == "Hinode"


def test_detector(sot_fg_wb_map):
    assert sot_fg_wb_map.detector == "SOT"


def test_reference_date(sot_fg_wb_map):
    assert sot_fg_wb_map.reference_date.isot == "1999-12-31T23:59:59.006"


def test_date(sot_fg_wb_map):
    assert sot_fg_wb_map.date.isot == "1999-12-31T23:59:59.006"


def test_private_date_setters(sot_fg_wb_map):
    _test_private_date_setters(sot_fg_wb_map)


def test_measurement(sot_fg_wb_map):
    assert sot_fg_wb_map.measurement is None


def test_instruments(sot_fg_wb_map):
    assert sot_fg_wb_map.Instruments == ['SOT/WB', 'SOT/NB', 'SOT/SP', 'SOT/CT']


def test_waves(sot_fg_wb_map):
    assert sot_fg_wb_map.Waves == ['6302A', 'BFI no move',
                                    'CN bandhead 3883', 'Ca II H line',
                                    'G band 4305', 'NFI no move', 'TF Fe I 6302',
                                    'TF Mg I 5172', 'TF Na I 5896',
                                    'blue cont 4504', 'green cont 5550',
                                    'red cont 6684']


def test_obstype(sot_fg_wb_map):
    assert sot_fg_wb_map.Observation_Type == ['FG (simple)',
                                              'FG focus scan', 'FG shuttered I and V',
                                              'FG shutterless I and V', 'FG shutterless I and V with 0.2s intervals',
                                              'FG shutterless Stokes', 'SP IQUV 4D array']


def test_wcs(sot_fg_wb_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='assuming Earth-based observer'):
        sot_fg_wb_map.pixel_to_world(0*u.pix, 0*u.pix)
