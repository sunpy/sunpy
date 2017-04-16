# -*- coding: utf-8 -*-
# Author: Rajul Srivastava<rajul@gmail.com>

#pylint: disable=W0613

from __future__ import absolute_import

from matplotlib import pyplot

from sunpy.visualization import visualization


__authors__ = ["Rajul Srivastava",]
__email__ = ["rajul09@gmail.com",]


def test_toggle_pylab():
	def dummy_function(dummy_arg='foo'):
		return dummy_arg

	assert pyplot.isinteractive() == False

	fn = visualization.toggle_pylab(dummy_function)
	assert fn == dummy_function
	assert fn() == 'foo'
	assert fn('bar') == 'bar'

	pyplot.ion()
	assert pyplot.isinteractive() == True
	
	fn_toggle = visualization.toggle_pylab(dummy_function)
	assert fn == dummy_function
	assert fn() == 'foo'
	assert fn('bar') == 'bar'
	assert pyplot.isinteractive() == True
