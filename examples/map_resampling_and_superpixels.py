# -*- coding: utf-8 -*-
"""
=========================================
Map Resampling and Superpixels
=========================================

In this example you will see how to resample a map using the resample method
(which implements interpolation) and superpixels.
"""

##############################################################################
# Start by importing the necessary modules.

from __future__ import print_function, division

import astropy.units as u

import matplotlib.pyplot as plt

import sunpy.map
import sunpy.data.sample

##############################################################################
# Sunpy sample data contains a number of suitable maps, where the sunpy.data.sample.NAME
# returns the location of the given FITS file.
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

##############################################################################
# This has a resolution of:
print(aia_map.dimensions)

##############################################################################
# To find out more specifics about this map and the instrument used, check it's
# metatdata:
print(aia_map.meta)

##############################################################################
# To reduce the angular resolution of the map you can use the `~sunpy.map.GenericMap.resample` method,
# specifying the dimensions as an Astropy Quantity in pixels:
dimensions = u.Quantity([40, 40], u.pixel)
aia_resampled_map = aia_map.resample(dimensions)
aia_resampled_map.peek(draw_limb=True, draw_grid=True)
# Note that this uses linear interpolation, you can change this using the method
# (‘neighbor’, ‘nearest’, ‘linear’ or ‘spline’) keyword argument option.

##############################################################################
# Similar to resampling you can use the `~sunpy.map.GenericMap.superpixel` method, this will reduce the
# resolution of the image by combining the number of pixels (in each dimension)
# in the dimensions argument into one single pixel.
# This can be used to increase the signal to noise ratio.
# For this the new dimensions must divide original image size exactly, for
# example you can reduce the AIA map resolution by a factor of 16 using:
dimensions = u.Quantity(aia_map.dimensions) / 16
aia_superpixel_map = aia_map.superpixel(dimensions)
aia_superpixel_map.peek(draw_limb=True, draw_grid=True)
plt.show()
