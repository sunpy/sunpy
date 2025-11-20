"""
Test cases for ASO-S HXIMap subclass.
"""
import pytest

import astropy.units as u
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.asos import HXIMap
from .helpers import _test_private_date_setters


@pytest.fixture(scope="module", params=[0,1])
def hxi_map(request):
    return get_dummy_map_from_header(get_test_filepath(f"hxi_imgcube_01e02t_20230501_130758_HXI_CLEAN_{request.param}.header"))


def test_HXIMap(hxi_map):
    """Tests the creation of HXIMap"""
    for tmp_map in hxi_map:
        assert isinstance(tmp_map, HXIMap)


def test_is_datasource_for(hxi_map):
    """Tests the is_datasource_for method of HXIMap."""
    for tmp_map in hxi_map:
        assert tmp_map.is_datasource_for(tmp_map.data, tmp_map.meta)


def test_observatory(hxi_map):
    """Tests the observatory property of the HXIMap object."""
    for tmp_map in hxi_map:
        assert tmp_map.observatory == "ASO-S"


def test_reference_date(hxi_map):
    for nth, tmp_map in enumerate(hxi_map):
        assert tmp_map.reference_date.isot == ("2023-05-01T13:07:58.713", '2023-05-01T13:08:16.962')[nth]


def test_date(hxi_map):
    for nth, tmp_map in enumerate(hxi_map):
        assert tmp_map.date.isot == ("2023-05-01T13:07:58.713", '2023-05-01T13:08:16.962')[nth]


def test_private_date_setters(hxi_map):
    for tmp_map in hxi_map:
        _test_private_date_setters(tmp_map)


def test_measurement(hxi_map):
    """Tests the measurement property of the HXIMap object."""
    for nth, tmp_map in enumerate(hxi_map):
        assert all(tmp_map.measurement == ([20, 30], [20, 30])[nth] * u.keV)


def test_norm_clip(hxi_map):
    # Tests that the default normalizer has clipping disabled
    for tmp_map in hxi_map:
        assert not tmp_map.plot_settings['norm'].clip


def test_new_instance_preserves_plot_settings(hxi_map):
    # Tests that the _new_instance method preserves the plot_settings
    # of the old instance. This is done on the HXI source as the HXIMap
    # constructor explicitly sets the cmap and norm and we want to test
    # that _new_instance persists the old custom plot_settings
    for tmp_map in hxi_map:
        tmp_map.plot_settings['norm'] = ImageNormalize(vmin=0.1, vmax=42)
        tmp_map.plot_settings['cmap'] = 'inferno'
        new_tmp_map = tmp_map._new_instance(tmp_map.data,
                                                  tmp_map.meta,
                                                  plot_settings=tmp_map.plot_settings)
        assert new_tmp_map.plot_settings['norm'].vmin == tmp_map.plot_settings['norm'].vmin
        assert new_tmp_map.plot_settings['norm'].vmax == tmp_map.plot_settings['norm'].vmax
        assert new_tmp_map.plot_settings['cmap'] == tmp_map.plot_settings['cmap']
        # If no plot settings are explicitly passed, the plot_settings should fall back to those
        # in the constructor
        new_tmp_map = tmp_map._new_instance(tmp_map.data, tmp_map.meta)
        assert new_tmp_map.plot_settings['norm'].vmin is None
        assert new_tmp_map.plot_settings['norm'].vmax is None
        assert new_tmp_map.plot_settings['cmap'] == new_tmp_map._get_cmap_name()


def test_wcs(hxi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
        for tmp_map in hxi_map:
            tmp_map.pixel_to_world(0*u.pix, 0*u.pix)
