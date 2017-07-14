# -*- coding: utf-8 -*-
# Author: Rajul Srivastava<rajul@gmail.com>

from __future__ import absolute_import

from matplotlib import pyplot

from sunpy.visualization import visualization


def test_toggle_pylab():
    def dummy_function(dummy_arg='foo'):
        return dummy_arg

    assert pyplot.isinteractive() is False

    fn = visualization.toggle_pylab(dummy_function)
    assert fn == dummy_function
    assert fn() == 'foo'
    assert fn('bar') == 'bar'

    pyplot.ion()
    assert pyplot.isinteractive() is True
