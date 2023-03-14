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

import matplotlib

###############################################################################
# We need to tell ``matplotlib`` to use a backend from ``mplcairo``.  The
# backend formally needs to be set prior to importing ``matplotlib.pyplot``.

print(f"Backend was: {matplotlib.get_backend()}")  # noqa
matplotlib.use("module://mplcairo.base")  # noqa
print(f"Backend is: {matplotlib.get_backend()}")  # noqa

###############################################################################
# We can now import everything else

print(f"Backend was: {matplotlib.get_backend()}")  # noqa
matplotlib.use("module://mplcairo.base")  # noqa
print(f"Backend is: {matplotlib.get_backend()}")  # noqa

import matplotlib.pyplot as plt

print(f"Backend was: {matplotlib.get_backend()}")  # noqa
matplotlib.use("module://mplcairo.base")  # noqa
print(f"Backend is: {matplotlib.get_backend()}")  # noqa

from mplcairo import operator_t

print(f"Backend is: {matplotlib.get_backend()}")  # noqa

import astropy.units as u

import sunpy.data.sample
import sunpy.map
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
# Plotting both the images, The clip_interval argument in both
# plot() functions sets the range of pixel values to display in the
# plot, with values outside of this range being clipped.

print(f"Backend was: {matplotlib.get_backend()}")
matplotlib.use("module://mplcairo.base")
print(f"Backend is: {matplotlib.get_backend()}")

fig = plt.figure()

print(f"Backend was: {matplotlib.get_backend()}")
matplotlib.use("module://mplcairo.base")
print(f"Backend is: {matplotlib.get_backend()}")

ax = fig.add_subplot(projection=a171)
_ = a171.plot(clip_interval=(1, 99.9995)*u.percent)
im131 = a131.plot(clip_interval=(1, 99.95)*u.percent)

# sphinx_gallery_defer_figures

###############################################################################
#  Now, patching the the ``im131`` plot with a screen operator.

print(f"Backend was: {matplotlib.get_backend()}")
matplotlib.use("module://mplcairo.base")
print(f"Backend is: {matplotlib.get_backend()}")

operator_t.SCREEN.patch_artist(im131)

# sphinx_gallery_defer_figures

###############################################################################
# Plot the result we get

print(f"Backend was: {matplotlib.get_backend()}")
matplotlib.use("module://mplcairo.base")
print(f"Backend is: {matplotlib.get_backend()}")

ax.set_title('mplcairo composite using screen blending')
plt.show()