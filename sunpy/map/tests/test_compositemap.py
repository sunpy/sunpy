# -*- coding: utf-8 -*-
"""
Test Composite Map
"""
import os
import pytest

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

import sunpy
import sunpy.map
import sunpy.coordinates
import sunpy.data.test
from sunpy.tests.helpers import figure_test
from sunpy.map import compositemap

testpath = sunpy.data.test.rootdir


@pytest.fixture
def aia171_test_map():
    return sunpy.map.Map(os.path.join(testpath, 'aia_171_level1.fits'))


@pytest.fixture
def heliographic_test_map():
    return sunpy.map.Map(os.path.join(testpath, 'heliographic_phase_map.fits.gz'))


@pytest.fixture
def composite_test_map():
    return sunpy.map.Map(aia171_test_map, heliographic_test_map, composite=True)


@figure_test
def test_plot_composite_map(composite_test_map):
    composite_test_map.plot()


@figure_test
def test_plot_composite_map_linewidths(composite_test_map):
    composite_test_map.plot(linewidths=4)
