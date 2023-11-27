"""
===============
Resampling Maps
===============

How to resample a map using the resample method, which implements interpolation, or
using superpixels, which combines pixels.
"""
import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.data.sample
import sunpy.map

##############################################################################
# We start with the sample data.

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

##############################################################################
# To reduce the angular resolution of the map, you can use the
# :meth:`~sunpy.map.GenericMap.resample` method, specifying the new dimensions
# in pixels. By default, this method uses linear interpolation but this can be
# changed with the ``method`` argument ('nearest', 'linear' or 'spline').

new_dimensions = [40, 40] * u.pixel
aia_resampled_map = aia_map.resample(new_dimensions)

##############################################################################
# Let's plot the result.

fig = plt.figure()
ax = fig.add_subplot(projection=aia_resampled_map)
aia_resampled_map.plot(axes=ax)
plt.show()

##############################################################################
# Another way to reduce the angular resolution of the map is by using the
# :meth:`~sunpy.map.GenericMap.superpixel` method, which combines pixels.
# The superpixel dimensions do not need to be square, and the intensity of
# each superpixel defaults to the sum of the constituent pixels. For example,
# you can reduce the AIA map resolution by a factor of 16 by specifying 16x16
# superpixels.

superpixel_size = [16, 16] * u.pixel
aia_superpixel_map = aia_map.superpixel(superpixel_size)

##############################################################################
# Let's plot the result.

fig = plt.figure()
ax = fig.add_subplot(projection=aia_superpixel_map)
aia_superpixel_map.plot(axes=ax)
plt.show()
