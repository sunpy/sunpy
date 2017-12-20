#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 13:32:41 2017

@author: ajit
"""

from __future__ import absolute_import
import os
import matplotlib.pyplot as plt
import pytest

import sunpy.map
from sunpy.visualization import mapcubeanimator
from sunpy.data import test
@pytest.fixture
def map_animator():
    testpath = test.rootdir
    aia_file = os.path.join(testpath, "aia_171_level1.fits")
    mapcube = sunpy.map.Map(aia_file)
    cube = sunpy.map.MapCube(mapcube)
    return mapcubeanimator.MapCubeAnimator(cube)

#test for mapcubeanimator instance
def test_mapcubeanimator_instances(map_animator):

    assert isinstance(map_animator, mapcubeanimator.MapCubeAnimator)
    assert map_animator.interval is 200
    assert map_animator.annotate is True
    assert map_animator.slider_ranges == [[0, 1]]
    assert map_animator.timer == None
    assert map_animator.active_slider == 0
    assert map_animator.button_func == []
    assert map_animator.remove_obj == []

def test_updatefig(map_animator):

    t = map_animator.updatefig(0, map_animator.im, 10)
    assert t == None

def test_plot_start_image(map_animator):
    t = map_animator.plot_start_image(plt.axes([0, 0, 0, 0]))
    s = map_animator.mapcube[0].plot(axes=plt.axes([0, 0, 0, 0]))
    assert t.axes == s.axes
    
