import os

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
            testpath, 'mdi_fd_Ic_6h_01d.5871.0000_s.fits'), sequence=True)


def test_construct_map_sequence_animator(test_map_sequence):
    map_animator = MapSequenceAnimator(test_map_sequence)
    assert isinstance(map_animator, MapSequenceAnimator)


@figure_test
def test_map_sequence_animator_wcs_simple_plot(test_map_sequence):
    map_animator = MapSequenceAnimator(test_map_sequence)
    map_animator._annotate_plot(0)
    return map_animator.fig


def test_axes():
    map_animator = MapSequenceAnimator(test_map_sequence)
    assert isinstance(map_animator, WCSAxes)
    start_img = map_animator.plot_start_image(map_animator.axes)
    assert isinstance(start_img.axes, WCSAxes)


@figure_test
def test_map_sequence_animator_wcs_update_plot(test_map_sequence):
    map_animator = MapSequenceAnimator(test_map_sequence)
    map_animator.updatefig(1, map_animator.im, map_animator.sliders[0]._slider)
    return map_animator.fig


@figure_test
def test_map_sequence_animator_wcs_colorbar_buttons(test_map_sequence):
    bf = [lambda x: x] * 10
    bl = ['button'] * 10
    map_animator = MapSequenceAnimator(
        test_map_sequence,
        colorbar=True,
        button_func=bf,
        button_labels=bl)
    map_animator.updatefig(
        1,
        map_animator.im,
        map_animator.sliders[0]._slider)
    return map_animator.fig


@figure_test
def test_map_sequence_animator_wcs_colorbar_buttons_default_labels(
        test_map_sequence):
    bf = [lambda x: x] * 10
    map_animator = MapSequenceAnimator(
        test_map_sequence, colorbar=True, button_func=bf)
    map_animator.updatefig(
        1,
        map_animator.im,
        map_animator.sliders[0]._slider)
    return map_animator.fig
