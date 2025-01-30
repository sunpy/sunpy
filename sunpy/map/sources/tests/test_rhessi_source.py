import pytest

import astropy.units as u
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.data.test import get_test_filepath
from sunpy.map import Map
from sunpy.map.sources.rhessi import RHESSIMap
from sunpy.util.exceptions import SunpyMetadataWarning
from .helpers import _test_private_date_setters


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


def test_new_instance_preserves_plot_settings(rhessi_map):
    # Tests that the _new_instance method preserves the plot_settings
    # of the old instance. This is done on the RHESSI source as the RHESSIMap
    # constructor explicitly sets the cmap and norm and we want to test
    # that _new_instance persists the old custom plot_settings
    rhessi_map.plot_settings['norm'] = ImageNormalize(vmin=0.1, vmax=42)
    rhessi_map.plot_settings['cmap'] = 'inferno'
    new_rhessi_map = rhessi_map._new_instance(rhessi_map.data,
                                              rhessi_map.meta,
                                              plot_settings=rhessi_map.plot_settings)
    assert new_rhessi_map.plot_settings['norm'].vmin == rhessi_map.plot_settings['norm'].vmin
    assert new_rhessi_map.plot_settings['norm'].vmax == rhessi_map.plot_settings['norm'].vmax
    assert new_rhessi_map.plot_settings['cmap'] == rhessi_map.plot_settings['cmap']
    # If no plot settings are explicitly passed, the plot_settings should fall back to those
    # in the constructor
    new_rhessi_map = rhessi_map._new_instance(rhessi_map.data, rhessi_map.meta)
    assert new_rhessi_map.plot_settings['norm'].vmin is None
    assert new_rhessi_map.plot_settings['norm'].vmax is None
    assert new_rhessi_map.plot_settings['cmap'] == new_rhessi_map._get_cmap_name()


def test_wcs(rhessi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning,
                      match='Missing metadata for observer: assuming Earth-based observer.*'):
        rhessi_map.pixel_to_world(0*u.pix, 0*u.pix)
