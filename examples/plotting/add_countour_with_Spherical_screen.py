"""
==========================================================
Add contours to a map using a SphericalScreen
==========================================================

How to do coordinate transformations using a SphericalScreen.

In this example, we will plot the 171 AIA image and add contours to it using `draw_contours`.
We will also plot the 193 AIA image and add contours to it. The contours of  193 AIA will drawn
by transforming the coordinates from the 171 image to the 193 image using the SphericalScreen class.
By specifying the axes of the 173 image in `draw_contours`, the contours will be drawn on the 193 image from
the 173 image.
"""

import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.data.sample
import sunpy.map
from sunpy.coordinates.screens import SphericalScreen

##############################################################################
# For this example, we will start with the sample data. To create the maps,
# we need the 171 and 193 AIA images. Both of these data can be downloaded
# with ``Fido``.

aia193 = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE)
aia171 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

##############################################################################
# Create a figure with two subplots. The first subplot will have the 171 image
# with contours drawn on it.
fig= plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(1, 2, 1, projection=aia171)
aia171.plot(axes=ax1, clip_interval=(1,99.9)*u.percent)
aia171.draw_contours(levels=[1,5]*u.percent, axes=ax1, colors='C0')
ax1.set_title("171 Contours on 171")

##############################################################################
# The second subplot will have the 193 image drawn on it.
ax2 = fig.add_subplot(1,2,2,projection=aia193)
aia193.plot(axes=ax2)

##############################################################################
# We will now draw the contours on the 193 image. To do this, we will use the
# `SphericalScreen` class. We will transform the coordinates from the 171 image
# to the 193 image using the `SphericalScreen` class. By specifying the
# axes of the 173 image in `draw_contours`, the contours will be drawn on the
# 193 image from the 171 image.
with SphericalScreen(aia171.observer_coordinate):
    aia171.draw_contours(levels=[1,5]*u.percent, axes=ax2, colors='C2')
ax2.set_title("171 Contours on 193")

##############################################################################
# Display the plot.
plt.show()
