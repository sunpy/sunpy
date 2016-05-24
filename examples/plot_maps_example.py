# -*- coding: utf-8 -*-
"""
=========================================
Interacting with Data Using SunPy Maps
=========================================

In this example you will be learning how to create and modify SunPy Map objects.
"""

##############################################################################
# Start by importing the necessary modules.

from __future__ import print_function, division

import numpy as np
import astropy.units as u

import sunpy.map
import sunpy.data.sample

##############################################################################
# SunPy Maps store 2D data in a numpy array and additional data in a metadata
# dictionary giving information relating to the data and instrument.
# You can create a Map in a number of ways, including loading a FITS file or URL:
# ``mymap = sunpy.map.Map('file1.fits')``
# ``mymap = sunpy.map.Map(url_str)``
# Or using creating manually by using tuple with the data/header within:
data = np.random.rand(20,15)
header = {}
manual_map = sunpy.map.Map((data, header))

##############################################################################
# The data numpy array and metadata dictionary can easily be accessed:
print(manual_map.data)
print(manual_map.meta)
# In this case notice that the metadata has been populated by default with the
# naxis details that correspond to the array used for the data.

##############################################################################
# You can quickly plot a map using the peek method:
manual_map.peek()

##############################################################################
# SunPy Maps have a number of attributes that can be accessed easily, such as
# the x and y ranges:
print(manual_map.xrange)
print(manual_map.yrange)
# These return astropy Quantity objects.
# In general the attributes are populated using details in the metadata and in
# this case there is no centre pixel or pixel size information given so SunPy
# is defaulting to assuming each pixel is 1 arcsec.
# This is in Helioprojective tangent projection in both longitude and latitude:
print(manual_map.coordinate_system)

##############################################################################
# A real map example is given in the sample data, where the sunpy.data.sample.NAME
# returns the location of the given FITS file.
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
aia_map.peek(draw_limb=True)

##############################################################################
# This has comprehensive metadata:
print(aia_map.meta)

##############################################################################
# Which allows it to accurately specify ranges:
print(aia_map.xrange)
print(aia_map.yrange)

# And find out information about the observation device and date:
print(aia_map.date)
print(aia_map.observatory)
print(aia_map.detector)
print(aia_map.exposure_time)
print(aia_map.coordinate_system)
print(aia_map.measurement)

##############################################################################
# To see only a part of the image you create a submap, specifying ranges in
# AstroPy Quantities:
rangex = u.Quantity([aia_map.xrange[0], 0 * u.arcsec])
rangey = u.Quantity([aia_map.yrange[0], 0 * u.arcsec])
aia_submap = aia_map.submap(rangex, rangey)
aia_submap.peek(draw_limb=True)

##############################################################################
# Similarly, if you want to reduce the angular resolution of the map you can use
# the `~sunpy.map.GenericMap.resample` method, specifying the dimensions as an Astropy Quantity in pixels:
dimensions = u.Quantity([50, 50], u.pixel)
aia_resampled_map = aia_map.resample(dimensions)
aia_resampled_map.peek(draw_limb=True, draw_grid=True)

##############################################################################
# Similar to resampling you can use the `~sunpy.map.GenericMap.superpixel` method, this will reduce the
# resolution of the image by combining the number of pixels (in each dimension)
# in the dimensions argument into one single pixel.
# This can be used to increase the signal to noise ratio.
# For this the new dimensions must divide original image size exactly.
dimensions = u.Quantity(aia_map.dimensions) / 16
aia_superpixel_map = aia_map.superpixel(dimensions)
aia_superpixel_map.peek(draw_limb=True)

##############################################################################
# Maps can also be rotated using the `~sunpy.map.GenericMap.rotate` method:
aia_rotated_submap = aia_submap.rotate(angle = 10 * u.deg)
aia_rotated_submap.peek(draw_limb=True, draw_grid=True)
# Note: the data array is expanded so that none of the original data is lost
# through clipping.