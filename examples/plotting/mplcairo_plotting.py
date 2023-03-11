"""
==========================================
Blending Plots using mplcairo
==========================================

This example will go through how you can blend two plots from SunPy
using mplcairo. ``matplotlib`` offers solely alpha-based transparency
for superimposing an image onto another, which can be restrictive
when trying to create visually appealing composite images. In
contrast, ``mplcairo`` provides a wide range of blending operators
for image overlays.
"""

import matplotlib; matplotlib.use("module://mplcairo.qt")
import matplotlib.pyplot as plt
from mplcairo import operator_t

import astropy.units as u

import sunpy.map
import sunpy.data.sample
from sunpy.coordinates import Helioprojective

###############################################################################
# Let's load two solar images, ``a171`` and ``a131`` we want to blend together
# and then reprojects the ``a131`` image onto the same coordinate system as the ``a171`` image

a171 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
a131 = sunpy.map.Map(sunpy.data.sample.AIA_131_IMAGE)
with Helioprojective.assume_spherical_screen(a171.observer_coordinate):
    a131 = a131.reproject_to(a171.wcs)

###############################################################################
# Now, using matplotlib and mplcairo to blend the two solar images.
# Setting axes according to ``a171`` image.

fig = plt.figure()
ax = fig.add_subplot(projection=a171)

###############################################################################
# Plotting both the images, The clip_interval argument in both 
# plot() functions sets the range of pixel values to display in the 
# plot, with values outside of this range being clipped.

_ = a171.plot(clip_interval=(1, 99.9995)*u.percent)
im131 = a131.plot(clip_interval=(1, 99.95)*u.percent)

###############################################################################
#  Now, patching the the ``im131`` plot with a screen operator.

operator_t.SCREEN.patch_artist(im131)

###############################################################################
# Plot the result we get

ax.set_title('mplcairo composite using screen blending')

plt.show()
