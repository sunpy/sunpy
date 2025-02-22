"""
===========================================
Overlaying off-disk contours using a screen
===========================================

This example shows how to use a screen to overlay off-disk contours from one
map onto another map.

When overlaying contours from one data set onto another data set, there is
almost always some small difference in the coordinate frame of the two data
sets, namely in the observation time and the observer location. With the
slight shift in perspective, any off-disk contours are not plotted by default
due to the lack of knowledge of where the contour exists in 3D space. (In
contrast, on-disk contours do not have this issue because the default
assumption is that on-disk contours lie on the surface of the Sun.) To change
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
# We use two AIA images that were not taken at exactly the same time. We
# crop the maps to an active region at the eastern limb.

aia171_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
aia193_map = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE)
aia171_map = aia171_map.submap((0, 500)*u.pix, top_right=(250, 750)*u.pix)
aia193_map = aia193_map.submap((0, 500)*u.pix, top_right=(250, 750)*u.pix)

##############################################################################
# First, let's plot the AIA 171 contours on their native map. The contours
# extend off disk due to the coronal loops.

fig = plt.figure()
ax = fig.add_subplot(projection=aia171_map)
aia171_map.plot(axes=ax)
aia171_map.draw_contours(axes=ax, levels=[1000, 2500] * u.DN, colors="green")

##############################################################################
# Next, let's try overlaying those same AIA 171 contours on the AIA 193 map.
# Note that none of the off-disk contours appear. Again, this is because it
# is not well defined where to plot such contours from the vantage point of
# the AIA 193 map.

fig = plt.figure()
ax = fig.add_subplot(projection=aia193_map)
aia193_map.plot(axes=ax, title="Default overlay behavior")
aia171_map.draw_contours(axes=ax, levels=[1000, 2500] * u.DN, colors="green")

##############################################################################
# Finally, let's overlay those contours again, but this time using the
# `~sunpy.coordinates.SphericalScreen` context manager. This screen assumes
# that off-disk points lie on the inside of a spherical screen. Note that the
# off-disk contours are now visible.

fig = plt.figure()
ax = fig.add_subplot(projection=aia193_map)
aia193_map.plot(axes=ax, title="Using SphericalScreen")
with SphericalScreen(aia171_map.observer_coordinate, only_off_disk=True):
    aia171_map.draw_contours(axes=ax, levels=[1000, 2500] * u.DN, colors="green")

plt.show()

# sphinx_gallery_thumbnail_number = 3
