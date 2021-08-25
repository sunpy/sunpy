import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.trace import TRACEMap
from sunpy.util.exceptions import SunpyMetadataWarning


@pytest.fixture(scope="module")
def trace_map():
    return get_dummy_map_from_header(get_test_filepath("tsi20010130_025823_a2.header"))


def test_fitstoTRACE(trace_map):
    """Tests the creation of TRACEMap using FITS."""
    assert isinstance(trace_map, TRACEMap)


def test_is_datasource_for(trace_map):
    """Test the is_datasource_for method of TRACEMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert trace_map.is_datasource_for(trace_map.data, trace_map.meta)


def test_measurement(trace_map):
    """Tests the measurement property of the TRACEMap object."""
    assert int(trace_map.measurement) == 171


def test_observatory(trace_map):
    """Tests the observatory property of the TRACEMap object."""
    assert trace_map.observatory == "TRACE"


def test_norm_clip(trace_map):
    # Tests that the default normalizer has clipping disabled
    assert not trace_map.plot_settings['norm'].clip


def test_wcs(trace_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        trace_map.pixel_to_world(0*u.pix, 0*u.pix)
