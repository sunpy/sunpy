"""
==============
Rotating a Map
==============

How to rotate a map.
"""
import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.data.sample
import sunpy.map

###############################################################################
# We start with the sample data
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

##############################################################################
# `~sunpy.map.GenericMap` provides the `~sunpy.map.GenericMap.rotate` method
# which accepts an angle. This returns a rotated map and does not rotate in
# place. The data array size is expanded so that none of the original data is
# lost due to clipping. Note that subsequent rotations are not compounded.
# The map is only rotated by the specified amount from the original maps
# orientation.
aia_rotated = aia_map.rotate(angle=30 * u.deg)

###############################################################################
# Let's now plot the results.
fig = plt.figure()
ax = plt.subplot(projection=aia_rotated)
aia_rotated.plot(clip_interval=(1, 99.99)*u.percent)
aia_rotated.draw_limb()
aia_rotated.draw_grid()
plt.show()
