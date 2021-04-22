"""Tests for EUI Solar Orbiter Map"""
import numpy as np
import pytest

import astropy.units as u
from astropy.io import fits

from sunpy.map import Map
from sunpy.map.sources import EUIMap


@pytest.fixture
def eui_map(scope="module"):
    raw_header = """
    """
    header = fits.Header.fromstring(raw_header, sep='\n')
    data = np.random.rand(header['naxis1'], header['naxis2'])
    return Map(data, header)


def test_EUIMap(eui_map):
    assert isinstance(eui_map, EUIMap)


def test_is_datasource_for(eui_map):
    assert eui_map.is_datasource_for(eui_map.data, eui_map.meta)


def test_observatory(eui_map):
    assert eui_map.observatory == "Solar Orbiter"


def test_measurement(eui_map):
    assert eui_map.measurement == u.Quantity(174, 'angstrom')


def test_exposure_time(eui_map):
    assert eui_map.exposure_time == u.Quantity(3, 's')


def test_level_number(eui_map):
    assert eui_map.processing_level == 2


def test_unit(eui_map):
    assert eui_map.unit == u.Unit('ct / s')


def test_norm_clip(eui_map):
    # Tests that the default normalizer has clipping disabled
    assert not eui_map.plot_settings['norm'].clip
