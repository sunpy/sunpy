"""
Test cases for Yohkoh SXTMap subclass
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.yohkoh import SXTMap
from sunpy.util.exceptions import SunpyMetadataWarning
from .helpers import _test_private_date_setters


@pytest.fixture
def sxt_map():
    with pytest.warns(SunpyMetadataWarning, match='Missing CTYPE'):
        return get_dummy_map_from_header(get_test_filepath("YohkohSXT.header"))


def test_fits_to_sxt(sxt_map):
    """Tests the creation of SXTMap using FITS."""
    assert isinstance(sxt_map, SXTMap)


def test_is_datasource_for(sxt_map):
    """Test the is_datasource_for method of SXTMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert sxt_map.is_datasource_for(sxt_map.data, sxt_map.meta)


def test_reference_date(sxt_map):
    assert sxt_map.reference_date.isot == "1991-11-05T11:10:24.018"


def test_date(sxt_map):
    assert sxt_map.date.isot == "1991-11-05T11:10:24.018"


def test_private_date_setters(sxt_map):
    _test_private_date_setters(sxt_map)


def test_observatory(sxt_map):
    """Tests the observatory property of the SXTMap object."""
    assert sxt_map.observatory == "Yohkoh"


def test_detector(sxt_map):
    """Tests the detector property of the SXTMap object"""
    assert sxt_map.detector == "SXT"


def test_measurement(sxt_map):
    """Tests the measurement property of the SXTMap object."""
    assert sxt_map.measurement == 'Al01'


def test_wavelength(sxt_map):
    """Tests that the wavelength of the SXTMap is always None"""
    assert sxt_map.wavelength is None


def test_heliographic_longitude(sxt_map):
    assert u.allclose(sxt_map.heliographic_longitude, 0.002357818089886795 * u.deg)


def test_heliographic_latitude(sxt_map):
    assert u.allclose(sxt_map.heliographic_latitude, 3.955251457847305 * u.deg)


def test_dsun(sxt_map):
    assert u.allclose(sxt_map.dsun, 148328341.63429108 * u.km)


def test_wcs(sxt_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing CTYPE'):
        sxt_map.pixel_to_world(0*u.pix, 0*u.pix)
