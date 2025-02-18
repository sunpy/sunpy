"""
===========================================
Overlaying off-disk contours using a screen
===========================================

This example shows how to use a screen to overlay off-disk contours from one
map onto another map.

When overlaying contours from one data set onto another data set, there is
almost always some small difference in the coordinate frame of the two data
sets, namely in the observation time and the observer location.  With the
slight shift in perspective, any off-disk contours are not plotted by default
due to the lack of knowledge of where the contour exists in 3D space.  (In
contrast, on-disk contours do not have this issue because the default
assumption is that on-disk contours lie on the surface of the Sun.)  To change
the behavior for off-disk contours, we can use a screen to specify an
assumption for the 3D location of 2D off-disk points.

See :ref:`sphx_glr_generated_gallery_map_transformations_reprojection_spherical_screen.py`
for a different use case for screens.
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
# We want to show how we can use `~sunpy.coordinates.SphericalScreen` as a
# way to transform the observer coordinates from one map to another. For this
# lets try drawing contours on the AIA 171 map(Plot 1) and use the observer
# coordinates from this for contours on the AIA 193 map(Plot 2). To do this,
# we will use the `~sunpy.coordinates.SphericalScreen` class so that it uses
# the observer coordinates of the AIA 171 image. So essentially we are
# transforming the coordinates of the contours for the AIA 193 image.
# If we dont use the `~sunpy.coordinates.SphericalScreen` class the contours
# are drawn on the 193 Map using the 193 Map's observer coordinates(Plot 3).

fig = plt.figure(figsize=(10, 4))

ax1 = fig.add_subplot(131, projection=aia171_map)
ax2 = fig.add_subplot(132, projection=aia193_map)
ax3 = fig.add_subplot(133, projection=aia193_map)

aia171_map.plot(axes=ax1, clip_interval=(1, 99.9) * u.percent)
aia171_map.draw_contours(levels=[1, 5] * u.percent, axes=ax1, colors='C0')
ax1.set_title("Contours on AIA 171 map")
ax1.set_xlabel("")

aia193_map.plot(axes=ax2)
with SphericalScreen(aia171_map.observer_coordinate):
    aia171_map.draw_contours(levels=[1, 5] * u.percent, axes=ax2, colors='C2')
ax2.set_title("Contours with SphericalScreen")
ax2.set_xlabel("")

aia193_map.plot(axes=ax3)
aia171_map.draw_contours(levels=[1, 5] * u.percent, axes=ax3, colors='C2')
ax3.set_title("Contours without SphericalScreen")
ax3.set_xlabel("")

for ax in [ax1, ax2, ax3]:
    ax.coords[0].set_ticks_visible(False)
    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[1].set_ticks_visible(False)
    ax.coords[1].set_ticklabel_visible(False)


plt.subplots_adjust(wspace=0)
fig.tight_layout()

plt.show()
