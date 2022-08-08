
import numpy as np
import pytest

from astropy.visualization.wcsaxes import WCSAxes

import sunpy.map
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test
from sunpy.visualization.animator import MapSequenceAnimator


@pytest.fixture
def test_map_sequence(aia171_test_map):
    return sunpy.map.Map(
        aia171_test_map,
        get_test_filepath('EIT/efz20040301.000010_s.fits'),
        sequence=True,
    )


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
    assert np.any(map1.data != map2.data)
