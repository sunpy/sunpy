#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 13:32:41 2017

@author: ajit
"""

from __future__ import absolute_import

import matplotlib.pyplot as plt
import pytest

import sunpy.map
from sunpy.visualization import mapcubeanimator
from sunpy.data import test


class map_animator(mapcubeanimator.MapCubeAnimator):
    pass

#test for mapcubeanimator instance
def test_mapcubeanimator_instances():
    file = test.get_test_filepath("FGMG4_20110214_030443.7.fits")
    dummy = sunpy.map.Map(file, cube=True)
    t = map_animator(dummy)
    assert isinstance(t, mapcubeanimator.MapCubeAnimator)

def test_get_main_axes():
    pass

def update_fig():
    file = test.get_test_filepath("FGMG4_20110214_030443.7.fits")
    dummy = sunpy.map.Map(file, cube=True)
    t = map_animator(dummy).updatefig(1)

    assert isinstance(t, mapcubeanimator.MapCubeAnimator.updatefig)
    assert t.val == 1
    assert t.im is None
    assert t.slider is None

def plot_image():
    file = test.get_test_filepath("FGMG4_20110214_030443.7.fits")
    dummy = sunpy.map.Map(file, cube=True)
    a = plt.axes([1, 2])
    t = map_animator(dummy).plot_start_image(a)

    assert isinstance(t, mapcubeanimator.MapCubeAnimator.plot_start_image)
    assert t.ax == a
