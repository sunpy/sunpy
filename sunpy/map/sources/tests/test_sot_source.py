import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.mapbase import SpatialPair
from sunpy.map.sources.hinode import SOTMap
from sunpy.util.exceptions import SunpyMetadataWarning
from .helpers import _test_private_date_setters

SOT_LOADABLE_HEADERS = {
    "HinodeSOT/HinodeSOT_FG_WB_G_band_4305_L0.header": {
        "instrument": "SOT/WB",
        "wave": "G band 4305",
        "date": "1999-12-31T23:59:59.006",
    },
    "HinodeSOT/HinodeSOT_FG_WB_CaIIH_L1.header": {
        "instrument": "SOT/WB",
        "wave": "Ca II H line",
        "date": "2006-11-06T12:00:32.887",
    },
    "HinodeSOT/HinodeSOT_FG_NB_TF_Na_I_5896_L1.header": {
        "instrument": "SOT/NB",
        "wave": "TF Na I 5896",
        "date": "2011-02-15T06:22:34.060",
    },
}


@pytest.fixture(scope="module", params=list(SOT_LOADABLE_HEADERS))
def sot_header(request):
    return request.param


@pytest.fixture(scope="module")
def sot_expected(sot_header):
    return SOT_LOADABLE_HEADERS[sot_header]


@pytest.fixture(scope="module")
def sot_map(sot_header):
    return get_dummy_map_from_header(get_test_filepath(sot_header))


def test_fitstoSOT(sot_map):
    assert isinstance(sot_map, SOTMap)


def test_is_datasource_for(sot_map):
    assert sot_map.is_datasource_for(sot_map.data, sot_map.meta)


def test_observatory(sot_map):
    assert sot_map.observatory == "Hinode"


def test_detector(sot_map):
    assert sot_map.detector == "SOT"


def test_instrument(sot_map, sot_expected):
    assert sot_map.instrument == sot_expected["instrument"]


def test_measurement(sot_map):
    assert sot_map.measurement is None


def test_wave_in_waves(sot_map, sot_expected):
    assert sot_expected["wave"] in sot_map._waves


def test_coordinate_system(sot_map):
    assert sot_map.coordinate_system == SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')


def test_reference_date(sot_map, sot_expected):
    assert sot_map.reference_date.isot == sot_expected["date"]


def test_date(sot_map, sot_expected):
    assert sot_map.date.isot == sot_expected["date"]


def test_private_date_setters(sot_map):
    _test_private_date_setters(sot_map)


def test_wcs(sot_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='assuming Earth-based observer'):
        sot_map.pixel_to_world(0*u.pix, 0*u.pix)
