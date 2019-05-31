# -*- coding: utf-8 -*-
"""
==============
Cropping a Map
==============

How to crop a map by using submap.
"""
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import sunpy.data.sample
import matplotlib.pyplot as plt

###############################################################################
# We start with the sample data
swap_map = sunpy.map.Map(sunpy.data.sample.SWAP_LEVEL1_IMAGE)

##############################################################################
# To crop the data you create a submap, specifying the top right and bottom
# left as SkyCoord objects.
top_right = SkyCoord(0 * u.arcsec, -200 * u.arcsec, frame=swap_map.coordinate_frame)
bottom_left = SkyCoord(-900 * u.arcsec, -900 * u.arcsec, frame=swap_map.coordinate_frame)
swap_submap = swap_map.submap(bottom_left, top_right)

###############################################################################
# Let's plot the results
fig = plt.figure()
ax = fig.add_subplot(111, projection=swap_submap)
swap_submap.plot()
swap_submap.draw_limb()
swap_submap.draw_grid()
plt.show()
