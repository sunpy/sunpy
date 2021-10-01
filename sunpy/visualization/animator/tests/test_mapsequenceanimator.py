import os

import numpy as np
import pytest

from astropy.visualization.wcsaxes import WCSAxes

import sunpy.data.test
import sunpy.map
from sunpy.tests.helpers import figure_test
from sunpy.visualization.animator import MapSequenceAnimator

testpath = sunpy.data.test.rootdir


@pytest.fixture
def test_map_sequence():
    return sunpy.map.Map(
        os.path.join(
            testpath, 'aia_171_level1.fits'), os.path.join(
            testpath, 'mdi.fd_Ic.20101015_230100_TAI.data.fits'), sequence=True)


def test_construct_map_sequence_animator(test_map_sequence):
    map_animator = MapSequenceAnimator(test_map_sequence)
    assert isinstance(map_animator, MapSequenceAnimator)


@figure_test
def test_map_sequence_animator_wcs_simple_plot(test_map_sequence):
    map_animator = MapSequenceAnimator(test_map_sequence)
    return map_animator.fig


def test_axes(test_map_sequence):
    map_animator = MapSequenceAnimator(test_map_sequence)
    assert isinstance(map_animator.axes, WCSAxes)


def test_map_sequence_animator_wcs_update_plot(test_map_sequence):
    map_animator = MapSequenceAnimator(test_map_sequence)
    map1 = map_animator.im.get_array()
    map_animator.updatefig(1, map_animator.im, map_animator.sliders[0]._slider)
    map2 = map_animator.im.get_array()
    assert np.all(map1.data != map2.data)
