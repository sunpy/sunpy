"""
Test cases for AIAMap subclass.
"""
import copy

import pytest

import astropy.units as u
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.data.test import get_test_filepath
from sunpy.map import Map
from sunpy.map.sources.sdo import AIAMap
from sunpy.tests.helpers import SKIP_GLYMUR
from .helpers import _test_private_date_setters

params = [get_test_filepath("aia_171_level1.fits")]
if not SKIP_GLYMUR:
    params += [get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2")]

__author__ = "Pritish C. (VaticanCameos)"


@pytest.fixture(scope="module", params=params)
def aia_map(request):
    return Map(request.param)


def test_aia_map(aia_map):
    """Tests the creation of AIAMap"""
    assert isinstance(aia_map, AIAMap)


def test_reference_date(aia_map):
    assert aia_map.reference_date.isot in ["2011-02-15T00:00:01.340", "2013-06-24T17:31:31.840"]


def test_date(aia_map):
    assert aia_map.date.isot in ["2011-02-15T00:00:00.340", "2013-06-24T17:31:30.840"]


def test_private_date_setters(aia_map):
    _test_private_date_setters(aia_map)


def test_is_datasource_for(aia_map):
    """Tests the is_datasource_for method of AIAMap."""
    assert aia_map.is_datasource_for(aia_map.data, aia_map.meta)


def test_observatory(aia_map):
    """Tests the observatory property of the AIAMap object."""
    assert aia_map.observatory == "SDO"


def test_measurement(aia_map):
    """Tests the measurement property of the AIAMap object."""
    # aiaimg has 171, jp2path has 193.
    assert aia_map.measurement.value in [171, 193]


def test_norm_clip(aia_map):
    # Tests that the default normalizer has clipping disabled
    assert not aia_map.plot_settings['norm'].clip


def test_new_instance_preserves_plot_settings(aia_map):
    # Tests that the _new_instance method preserves the plot_settings
    # of the old instance. This is done on the AIA source as the AIAMap
    # constructor explicitly sets the cmap and norm and we want to test
    # that _new_instance persists the old custom plot_settings
    aia_map.plot_settings['norm'] = ImageNormalize(vmin=0.1, vmax=42)
    aia_map.plot_settings['cmap'] = 'inferno'
    new_aia_map = aia_map._new_instance(aia_map.data,
                                        aia_map.meta,
                                        plot_settings=aia_map.plot_settings)
    assert new_aia_map.plot_settings['norm'].vmin == aia_map.plot_settings['norm'].vmin
    assert new_aia_map.plot_settings['norm'].vmax == aia_map.plot_settings['norm'].vmax
    assert new_aia_map.plot_settings['cmap'] == aia_map.plot_settings['cmap']
    # If no plot settings are explicitly passed, the plot_settings should fall back to those
    # in the constructor
    new_aia_map = aia_map._new_instance(aia_map.data, aia_map.meta)
    assert new_aia_map.plot_settings['norm'].vmin is None
    assert new_aia_map.plot_settings['norm'].vmax is None
    assert new_aia_map.plot_settings['cmap'] == new_aia_map._get_cmap_name()


def test_wcs(aia_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    aia_map.pixel_to_world(0*u.pix, 0*u.pix)


def test_missing_tobs(aia_map):
    # Should fall back to base reference_date for reference date if T_OBS is missing
    new_meta = copy.deepcopy(aia_map.meta)
    new_meta.pop('T_OBS')
    new_aia_map = aia_map._new_instance(aia_map.data, new_meta)
    assert new_aia_map.reference_date == super(type(new_aia_map), new_aia_map).reference_date
