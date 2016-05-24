# -*- coding: utf-8 -*-
"""
=========================================
Map Rotation
=========================================

In this example we rotate a map.
"""

##############################################################################
# Start by importing the necessary modules.

from __future__ import print_function, division

import astropy.units as u

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
# Maps can also be rotated by using the `~sunpy.map.GenericMap.rotate` method
# with a specified angle, supplied as an Astropy Quantity with angular units:
aia_rotated = aia_map.rotate(angle = 30 * u.deg)
aia_rotated.peek(draw_limb=True, draw_grid=True)
# Or using Radians
aia_rotated = aia_map.rotate(angle = 0.5 * u.rad)
# Note: the data array is expanded so that none of the original data is lost
# through clipping.
# Also note that subsequent rotations are not compunded, the map is only rotated
# by the specified amount from the original maps orientation.