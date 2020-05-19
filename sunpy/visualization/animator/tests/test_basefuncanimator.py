import warnings
from functools import partial
from itertools import product

import matplotlib.animation as mplanim
import matplotlib.axes as maxes
import matplotlib.backend_bases as mback
import matplotlib.pyplot as plt
import numpy as np
import pytest

import astropy.units as u
import astropy.wcs

import sunpy.data.test
import sunpy.map
from sunpy.tests.helpers import figure_test
from sunpy.time import parse_time
from sunpy.util.exceptions import SunpyDeprecationWarning
from sunpy.visualization.animator import ArrayAnimator, BaseFuncAnimator, ImageAnimatorWCS, LineAnimator, base


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
    ranges = [(0, 3), (0, 3)]

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

    event.key = 'up'
    tfa._key_press(event)
    assert tfa.active_slider == 1

    tfa.slider_buttons[tfa.active_slider]._button.clicked = False
    event.key = 'p'
    tfa._click_slider_button(event=event, button=tfa.slider_buttons[tfa.active_slider]._button,
                             slider=tfa.sliders[tfa.active_slider]._slider)
    assert tfa.slider_buttons[tfa.active_slider]._button.label._text == "||"

    tfa._key_press(event)
    assert tfa.slider_buttons[tfa.active_slider]._button.label._text == ">"

    event.key = 'left'
    tfa._key_press(event)
    assert tfa.sliders[tfa.active_slider]._slider.val == tfa.sliders[tfa.active_slider]._slider.valmax

    event.key = 'right'
    tfa._key_press(event)
    assert tfa.sliders[tfa.active_slider]._slider.val == tfa.sliders[tfa.active_slider]._slider.valmin

    event.key = 'right'
    tfa._key_press(event)
    assert tfa.sliders[tfa.active_slider]._slider.val == tfa.sliders[tfa.active_slider]._slider.valmin + 1

    event.key = 'left'
    tfa._key_press(event)
    assert tfa.sliders[tfa.active_slider]._slider.val == tfa.sliders[tfa.active_slider]._slider.valmin

    tfa._start_play(event, tfa.slider_buttons[tfa.active_slider]._button,
                    tfa.sliders[tfa.active_slider]._slider)
    assert tfa.timer

    tfa._stop_play(event)
    assert tfa.timer is None

    tfa._slider_changed(val=2, slider=tfa.sliders[tfa.active_slider]._slider)
    assert np.array(tfa.im.get_array()).all() == data[2].all()

    event.inaxes = tfa.sliders[0]
    tfa._mouse_click(event)
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
        super().update_plot(val, artist, slider)


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
    aanim = ArrayAnimatorTest(data=data)
    out_axis_ranges, out_extent = aanim._sanitize_axis_ranges(axis_ranges=axis_ranges,
                                                              data_shape=data_shape)
    assert exp_extent == out_extent
    assert np.array_equal(exp_axis_ranges[1], out_axis_ranges[1])
    assert callable(out_axis_ranges[0])
    assert np.array_equal(exp_axis_ranges[0], out_axis_ranges[0](np.arange(10)))


xdata = np.tile(np.linspace(0, 100, 11), (5, 5, 1))


@pytest.mark.parametrize('plot_axis_index, axis_ranges, xlabel, xlim',
                         [(-1, None, None, None),
                          (-1, [None, None, xdata], 'x-axis', None)])
def test_lineanimator_init(plot_axis_index, axis_ranges, xlabel, xlim):
    data = np.random.random((5, 5, 10))
    LineAnimator(data=data, plot_axis_index=plot_axis_index, axis_ranges=axis_ranges,
                 xlabel=xlabel, xlim=xlim)


@figure_test
def test_lineanimator_figure():
    np.random.seed(1)
    data_shape0 = (10, 20)
    data0 = np.random.rand(*data_shape0)
    plot_axis0 = 1
    slider_axis0 = 0
    xdata = np.tile(np.linspace(
        0, 100, (data_shape0[plot_axis0] + 1)), (data_shape0[slider_axis0], 1))

    ani = LineAnimator(data0, plot_axis_index=plot_axis0, axis_ranges=[None, xdata])

    return ani.fig


@figure_test
def test_imageanimator_figure():
    AIA_171 = sunpy.data.test.get_test_filepath('aia_171_level1.fits')
    KCOR = sunpy.data.test.get_test_filepath('20181209_180305_kcor_l1.5_rebinned.fits')
    map_seuence = sunpy.map.Map(AIA_171, KCOR, sequence=True)
    sequence_array = map_seuence.as_array()
    wcs_input_dict = {f'{key}{n+1}': map_seuence.all_meta()[0].get(f'{key}{n}')
                      for n, key in product([1, 2], ['CTYPE', 'CUNIT', 'CDELT'])}
    t0, t1 = map(parse_time, [k['date-obs'] for k in map_seuence.all_meta()])
    time_diff = (t1 - t0).to(u.s)
    wcs_input_dict.update(
        {'CTYPE1': 'Time', 'CUNIT1': time_diff.unit.name, 'CDELT1': time_diff.value})
    wcs = astropy.wcs.WCS(wcs_input_dict)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SunpyDeprecationWarning)
        wcs_anim = ImageAnimatorWCS(sequence_array, wcs=wcs, vmax=1000, image_axes=[0, 1])

    return wcs_anim.fig
