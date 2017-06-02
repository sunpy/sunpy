# -*- coding: utf-8 -*-
"""
=========================================
Submaps and Cropping
=========================================

In this example we demonstrate how to get a submap of a map.
"""

##############################################################################
# Start by importing the necessary modules.

from __future__ import print_function, division

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import sunpy.data.sample

##############################################################################
# Sunpy sample data contains a number of suitable maps, where the sunpy.data.sample.NAME
# returns the location of the given FITS file.
swap_map = sunpy.map.Map(sunpy.data.sample.SWAP_LEVEL1_IMAGE)

##############################################################################
# This has resolution and ranges of:
print(swap_map.dimensions)
print(swap_map.data)

##############################################################################
# To find out more specifics about this map and the instrument used, check it's
# metatdata:
print(swap_map.meta)

##############################################################################
# To crop the data you create a submap, specifying the top right and bottom
# left as SkyCoord objects.
top_right = SkyCoord(0*u.arcsec, -200 * u.arcsec, frame=swap_map.coordinate_frame)
bottom_left = SkyCoord(-900 * u.arcsec, -900 * u.arcsec, frame=swap_map.coordinate_frame)
swap_submap = swap_map.submap(bottom_left, top_right)
swap_submap.peek(draw_limb=True, draw_grid=True)
