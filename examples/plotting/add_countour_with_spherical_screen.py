"""
=============================================
Add contours to a map using a SphericalScreen
=============================================

In this example, we try to show that SphericalScreen can be used as a general method
to transform the coordinates(change the frame of reference i.e. the observer)
from one image to another as seen by a new observer, even without reprojection.
"""

import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.data.sample
import sunpy.map
from sunpy.coordinates.screens import SphericalScreen

##############################################################################
# For this example, we will start with the sunpy sample data.

aia193_map = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE)
aia171_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

##############################################################################
# Create a figure with two subplots. The first subplot will have the 171 image
# with contours drawn on it. The second subplot will have the 193 image drawn
# on it.
fig= plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(1, 3, 1, projection=aia171_map)
aia171_map.plot(axes=ax1, clip_interval=(1,99.9)*u.percent)
aia171_map.draw_contours(levels=[1,5]*u.percent, axes=ax1, colors='C0')
ax1.set_title("Contours on 173")
ax2 = fig.add_subplot(1,3,2,projection=aia193_map)
aia193_map.plot(axes=ax2)

##############################################################################
# We will now draw the contours on the AIA 193 map. To do this, we will use the
# `SphericalScreen` class so that it uses the observer/coordinates of the AIA 171
# image. So essentially we are transforming the coordinates of the contours
# for the AIA 193 image.

with SphericalScreen(aia171_map.observer_coordinate):
    aia171_map.draw_contours(levels=[1, 5]*u.percent, axes=ax2, colors='C2')
ax2.set_title("Contours on 193 with SphericalScreen")

##############################################################################
# without `SphereicalScreen` class the contours are drawn on the 193 image
# using the 193 image's observer/coordinates.

ax3 = fig.add_subplot(1,3,3,projection=aia193_map)
aia193_map.plot(axes=ax3)
aia171_map.draw_contours(levels=[1,5]*u.percent, axes=ax3, colors='C2')
ax3.set_title("Contours on 193 without SphericalScreen")

plt.show()
