# -*- coding: utf-8 -*-

from functools import partial
import pytest

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as mplanim
import matplotlib.backend_bases as mback

from sunpy.visualization.animator import base, BaseFuncAnimator


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
    event = mback.KeyEvent(name='key_press_event', canvas=fig.canvas, key='down')
    tfa._key_press(event)
    assert tfa.active_slider == 0


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


def test_edges_to_centers_nd():
    edges_axis = 0
    axis_range = np.zeros((10, 2))
    axis_range[:, 0] = np.arange(10, 20)
    expected = np.zeros((9,2))
    expected[:, edges_axis] = np.arange(10.5, 19)
    output = base.edges_to_centers_nd(axis_range, edges_axis)
    assert np.array_equal(output, expected)
