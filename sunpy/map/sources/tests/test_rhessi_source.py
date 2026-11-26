import pytest

import astropy.units as u

from sunpy.data.test import get_test_filepath
from sunpy.map import Map
from sunpy.map.sources.rhessi import RHESSIMap
from sunpy.util.exceptions import SunpyDeprecationWarning
from .helpers import _test_private_date_setters

# The observer is now resolved when the map is constructed rather than on first
# use of the WCS, so these maps emit the warning during fixture setup.
pytestmark = pytest.mark.filterwarnings("ignore:Missing metadata for observer")


@pytest.fixture(scope="module")
def rhessi_map():
    return Map(get_test_filepath("hsi_image_20101016_191218.fits"))


def test_RHESSIMap(rhessi_map):
    """Tests the creation of RHESSIMap"""
    assert isinstance(rhessi_map, RHESSIMap)


def test_is_datasource_for(rhessi_map):
    """Tests the is_datasource_for method of RHESSIMap."""
    assert rhessi_map.is_datasource_for(rhessi_map.data, rhessi_map.meta)


def test_observatory(rhessi_map):
    """Tests the observatory property of the RHESSIMap object."""
    assert rhessi_map.observatory == "RHESSI"


def test_reference_date(rhessi_map):
    assert rhessi_map.reference_date.isot == "2010-10-16T19:12:18.000"


def test_date(rhessi_map):
    assert rhessi_map.date.isot == "2010-10-16T19:12:18.000"


def test_private_date_setters(rhessi_map):
    _test_private_date_setters(rhessi_map)


def test_measurement(rhessi_map):
    """Tests the measurement property of the RHESSIMap object."""
    assert all(rhessi_map.measurement == [12, 25] * u.keV)


def test_norm_clip(rhessi_map):
    # Tests that the default normalizer has clipping disabled
    assert not rhessi_map.plot_settings['norm'].clip


def test_wcs(rhessi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    rhessi_map.wcs.pixel_to_world(0, 0)


def test_old_sig_sources(rhessi_map):
    with pytest.warns(SunpyDeprecationWarning,
                      match="header as a positional argument is deprecated"):
        smap = RHESSIMap(rhessi_map.data, rhessi_map.meta)
    assert isinstance(smap, RHESSIMap)
