# -*- coding: utf-8 -*-

from functools import partial
import pytest

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import widgets
import matplotlib.axes as maxes
import matplotlib.animation as mplanim
import matplotlib.backend_bases as mback

from sunpy.visualization.animator import base, BaseFuncAnimator, ArrayAnimator, LineAnimator


class FuncAnimatorTest(BaseFuncAnimator):
    def plot_start_image(self, ax):
        im = ax.imshow(self.data[0])
        if self.if_colorbar:
            self._add_colorbar(im)
        return im


def update_plotval(val, im, slider, data):
    i = int(val)
    im.set_array(data[i])


def button_func1(*args, **kwargs):
    print(*args, **kwargs)


@pytest.mark.parametrize('fig, colorbar, buttons', ((None, False, [[], []]),
                                                    (plt.figure(), True, [[button_func1], ["hi"]])))
def test_base_func_init(fig, colorbar, buttons):
    data = np.random.random((3, 10, 10))
    func0 = partial(update_plotval, data=data)
    func1 = partial(update_plotval, data=data*10)
    funcs = [func0, func1]
    ranges = [(0, 3), (4, 7)]

    tfa = FuncAnimatorTest(data, funcs, ranges, fig=fig, colorbar=colorbar,
                           button_func=buttons[0],
                           button_labels=buttons[1])

    tfa.label_slider(0, "hello")
    assert tfa.sliders[0]._slider.label.get_text() == "hello"

    tfa._set_active_slider(1)
    assert tfa.active_slider == 1

    fig = plt.figure()
    event = mback.KeyEvent(name='key_press_event', canvas=fig.canvas, key='down', valinit=0)
    tfa._key_press(event)
    assert tfa.active_slider == 0

    slide = tfa.active_slider
    slider = widgets.Slider(fig.gca(), str(slide), ranges[slide][0], ranges[slide][1])
    slider.slider_ind = slide
    butt = widgets.Button(fig.gca(), ">")
    butt.clicked = False
    event.key = 'p'
    tfa._click_slider_button(event=event, button=butt, slider=slider)
    assert butt.label._text == "||"

    tfa._start_play(event, butt, slider)
    assert tfa.timer
    
    tfa._stop_play(event)
    assert tfa.timer == None

    tfa._previous(slider)
    assert slider.val == slider.valmin

    tfa._step(slider)
    assert slider.val == slider.valmax



@pytest.fixture
def funcanimator():
    data = np.random.random((3, 10, 10))
    func = partial(update_plotval, data=data)
    funcs = [func]
    ranges = [(0, 3)]

    return FuncAnimatorTest(data, funcs, ranges)


def test_to_anim(funcanimator):
    ani = funcanimator.get_animation()
    assert isinstance(ani, mplanim.FuncAnimation)


def test_to_axes(funcanimator):
    ax = funcanimator._get_main_axes()
    assert isinstance(ax, maxes._subplots.SubplotBase)


def test_edges_to_centers_nd():
    edges_axis = 0
    axis_range = np.zeros((10, 2))
    axis_range[:, 0] = np.arange(10, 20)
    expected = np.zeros((9, 2))
    expected[:, edges_axis] = np.arange(10.5, 19)
    output = base.edges_to_centers_nd(axis_range, edges_axis)
    assert np.array_equal(output, expected)


class ArrayAnimatorTest(ArrayAnimator):
    def __init__(self, data):
        self.naxis = data.ndim
        self.image_axes = [1]
        self.slider_axes = [0]

    def plot_start_image(self, ax):
        pass

    def update_plot(self, val, artist, slider):
        pass


axis_ranges1 = np.tile(np.linspace(0, 100, 21), (10, 1))


@pytest.mark.parametrize('axis_ranges, exp_extent, exp_axis_ranges',
                          [([None, None], [-0.5, 19.5],
                          [np.arange(10), np.array([-0.5, 19.5])]),

                          ([[0, 10], [0, 20]], [0, 20],
                          [np.arange(0.5, 10.5), np.asarray([0, 20])]),

                          ([np.arange(0, 11), np.arange(0, 21)], [0, 20],
                          [np.arange(0.5, 10.5), np.arange(0.5, 20.5)]),
                          
                          ([None, axis_ranges1], [0.0, 100.0],
                          [np.arange(10), base.edges_to_centers_nd(axis_ranges1, 1)])])
def test_sanitize_axis_ranges(axis_ranges, exp_extent, exp_axis_ranges):
    data_shape = (10, 20)
    data = np.random.rand(*data_shape)
    edges_axis = 1
    aanim = ArrayAnimatorTest(data=data)
    out_axis_ranges, out_extent = aanim._sanitize_axis_ranges(axis_ranges=axis_ranges,
                                                                data_shape=data_shape)
    assert exp_extent == out_extent
    assert np.array_equal(exp_axis_ranges[0], out_axis_ranges[0])
    assert np.array_equal(exp_axis_ranges[1], out_axis_ranges[1])


xdata = np.tile(np.linspace(0, 100, 11), (5, 5, 1))


@pytest.mark.parametrize('plot_axis_index, axis_ranges, xlabel, xlim',
                        [(-1, None, None, None),
                        (-1, [None, None, xdata], 'x-axis', None)])
def test_lineanimator_init(plot_axis_index, axis_ranges, xlabel, xlim):
    data = np.random.random((5, 5, 10))
    ani = LineAnimator(data=data, plot_axis_index=plot_axis_index, axis_ranges=axis_ranges,
                      xlabel=xlabel, xlim=xlim)
