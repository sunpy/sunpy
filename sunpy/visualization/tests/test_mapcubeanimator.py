#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 13:32:41 2017

@author: ajit
"""

from __future__ import absolute_import

import matplotlib.pyplot as plt

import pytest

import sunpy

from sunpy.visualization import mapcubeanimator

@pytest.fixture
def mapcube_animator():
    mapcube = sunpy.map.Map()
    return mapcubeanimator.MapCubeAnimator(mapcube, slider_ranges = [[0, 5]])

#test for mapcubeanimator instances
def test_mapcubeanimator_instances(mapcube_animator):
    t = sunpy.map.Map()
    s = plt.axes([0.1, 0.8, 0.8, 0.1])
    assert mapcube_animator.mapcube == t
    assert mapcube_animator.annotate is True
    assert mapcube_animator.colorbar is True
    assert mapcube_animator.fig == s
    assert mapcube_animator.button_labels is None
    assert mapcube_animator.button_func is None
    assert mapcube_animator.interval == 1

def test_mapcubeanimator_plot_start_image(mapcube_animator):
    pass

def test_mapcubeanimator_updatefig():
    pass
