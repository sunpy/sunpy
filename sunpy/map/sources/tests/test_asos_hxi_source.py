"""
Test cases for ASO-S HXIMap subclass.
"""
import pytest

import astropy.units as u
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.asos import HXIMap
from sunpy.util.exceptions import SunpyMetadataWarning
from .helpers import _test_private_date_setters


@pytest.fixture(scope="module", params=[0,1])
def hxi_map(request):
    return get_dummy_map_from_header(get_test_filepath(f"hxi_imgcube_01e02t_20230501_130758_HXI_CLEAN_{request.param}.header"))


def test_HXIMap(hxi_map):
    """Tests the creation of HXIMap"""
    assert isinstance(hxi_map, HXIMap)


def test_is_datasource_for(hxi_map):
    """Tests the is_datasource_for method of HXIMap."""
    assert hxi_map.is_datasource_for(hxi_map.data, hxi_map.meta)


def test_observatory(hxi_map):
    """Tests the observatory property of the HXIMap object."""
    assert hxi_map.observatory == "ASO-S"


def test_reference_date(hxi_map):
    assert hxi_map.reference_date.isot in ["2023-05-01T13:07:58.713", "2023-05-01T13:08:16.962"]


def test_date(hxi_map):
    assert hxi_map.date.isot in ["2023-05-01T13:07:58.713", "2023-05-01T13:08:16.962"]


def test_unit(hxi_map):
    assert hxi_map.unit == u.photon / (u.s * u.arcsec**2 * u.cm**2)


def test_broken_units(hxi_map):
    # Check that no unit or bunit leads to no unit
    hxi_map.meta["unit"] = None
    assert hxi_map.unit is None

    # Check that other units are parsed
    hxi_map.meta["unit"] = "second"
    assert hxi_map.unit == u.s

    # Check that it falls back to fits standard key
    hxi_map.meta["bunit"] = "second"
    assert hxi_map.unit == u.s

    # Check that non-fits units still raise
    hxi_map.meta["unit"] = "spam"
    with pytest.warns(SunpyMetadataWarning):
        hxi_map.unit


def test_private_date_setters(hxi_map):
    _test_private_date_setters(hxi_map)


def test_measurement(hxi_map):
    """Tests the measurement property of the HXIMap object."""
    assert u.allclose(hxi_map.measurement, [20, 30] * u.keV)


def test_norm_clip(hxi_map):
    # Tests that the default normalizer has clipping disabled
    assert not hxi_map.plot_settings['norm'].clip


def test_new_instance_preserves_plot_settings(hxi_map):
    # Tests that the _new_instance method preserves the plot_settings
    # of the old instance. This is done on the HXI source as the HXIMap
    # constructor explicitly sets the cmap and norm and we want to test
    # that _new_instance persists the old custom plot_settings
    hxi_map.plot_settings['norm'] = ImageNormalize(vmin=0.1, vmax=42)
    hxi_map.plot_settings['cmap'] = 'inferno'
    new_hxi_map = hxi_map._new_instance(hxi_map.data,
                                        hxi_map.meta,
                                        plot_settings=hxi_map.plot_settings)
    assert new_hxi_map.plot_settings['norm'].vmin == hxi_map.plot_settings['norm'].vmin
    assert new_hxi_map.plot_settings['norm'].vmax == hxi_map.plot_settings['norm'].vmax
    assert new_hxi_map.plot_settings['cmap'] == hxi_map.plot_settings['cmap']
    # If no plot settings are explicitly passed, the plot_settings should fall back to those
    # in the constructor
    new_hxi_map = hxi_map._new_instance(hxi_map.data, hxi_map.meta)
    assert new_hxi_map.plot_settings['norm'].vmin is None
    assert new_hxi_map.plot_settings['norm'].vmax is None
    assert new_hxi_map.plot_settings['cmap'] == new_hxi_map._get_cmap_name()


def test_wcs(hxi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    hxi_map.pixel_to_world(0*u.pix, 0*u.pix)
