"""
Test cases for AIAMap subclass.
"""

import pytest

import astropy.units as u

from sunpy.data.test import get_test_filepath
from sunpy.map import Map
from sunpy.map.sources.sdo import AIAMap
from sunpy.tests.helpers import SKIP_GLYMUR

params = [get_test_filepath("aia_171_level1.fits")]
if not SKIP_GLYMUR:
    params += [get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2")]

__author__ = "Pritish C. (VaticanCameos)"


@pytest.fixture(scope="module", params=params)
def aia_map(request):
    return Map(request.param)


def test_AIAMap(aia_map):
    """Tests the creation of AIAMap"""
    assert isinstance(aia_map, AIAMap)


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


def test_wcs(aia_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    aia_map.pixel_to_world(0*u.pix, 0*u.pix)
