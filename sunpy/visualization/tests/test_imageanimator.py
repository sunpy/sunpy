# -*- coding: utf-8 -*-
# Author: Rajul Srivastava <rajul@gmail.com>

from __future__ import absolute_import

import matplotlib.pyplot as plt

import pytest

from sunpy.visualization import imageanimator


@pytest.fixture
def slider_pb():
    test_slider_ax = plt.axes([0.1, 0.8, 0.8, 0.1])
    return imageanimator.SliderPB(test_slider_ax, 'foo', -5, 5)


# Tests for the SliderPB class
def test_SliderPB_instance(slider_pb):
    dummy_slider_ax = plt.axes([0.1, 0.8, 0.8, 0.1])

    assert slider_pb.ax == dummy_slider_ax
    assert slider_pb.label.get_text() == 'foo'
    assert slider_pb.valmin == -5
    assert slider_pb.valmax == 5
    assert slider_pb.valinit == 0.5
    assert slider_pb.valfmt == '%1.2f'
    assert slider_pb.closedmin is True
    assert slider_pb.closedmax is True
    assert slider_pb.slidermin is None
    assert slider_pb.slidermax is None

    assert slider_pb.changed_args == {}

    assert slider_pb.cnt == 0


def test_SliderPB_set_val(slider_pb):
    pass


def test_SliderPB_on_changed(slider_pb):
    def dummy_function():
        pass

    t = slider_pb.on_changed(dummy_function)

    assert t == 0
    assert slider_pb.cnt == 1
    assert len(slider_pb.observers) == 1
    assert slider_pb.observers[t] == dummy_function
    assert len(slider_pb.changed_args) == 1
    assert slider_pb.changed_args[t] == ()
