"""Tests for EUI Solar Orbiter Map"""
import pytest

from sunpy.map import Map
from sunpy.map.sources import EUIMap


@pytest.fixture
def eui_map():
    return Map(eui_test_fits_file)


def test_EUIMap(eui_map):
    assert isinstance(eui_map, EUIMap)


def test_is_datasource_for(eui_map):
    assert eui_map.is_datasource_for(eui_map.data, eui_map.meta)


def test_observatory(eui_map):
    assert eui_map.observatory == "Solar Orbiter"


def test_measurement(eui_map):
    assert eui_map.measurement == u.Quantity(174, 'angstrom')


def test_norm_clip(eui_map):
    # Tests that the default normalizer has clipping disabled
    assert not eui_map.plot_settings['norm'].clip
