"""
================
Cropping a Map
================

How to crop a map to a set of world coordinates with
:meth:`~ndcube.NDCube.crop`.

:meth:`~ndcube.NDCube.crop` returns the smallest rectangle of pixels containing
all the coordinates it is given. Unlike :meth:`~sunpy.map.GenericMap.submap`, it
does not assume that the coordinates are opposite corners of a rectangle, and
it accepts any number of coordinates. For a comparison of the two methods, see
:ref:`sunpy-topic-guide-submap-vs-crop`.
"""
# sphinx_gallery_tags = ["Map", "AIA", "Coordinates"]

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map

###############################################################################
# We start with the sample data.

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

###############################################################################
# With two coordinates, :meth:`~ndcube.NDCube.crop` returns the smallest
# region of the map that contains both.

bottom_left = SkyCoord(-800 * u.arcsec, -300 * u.arcsec, frame=aia_map.coordinate_frame)
top_right = SkyCoord(-200 * u.arcsec, 400 * u.arcsec, frame=aia_map.coordinate_frame)
two_point_crop = aia_map.crop(bottom_left, top_right)

###############################################################################
# You can pass any number of coordinates. This is useful when you want to cover
# several features rather than cut out a box. Here we ask for the region
# containing five scattered points.

points = SkyCoord([-600, -100, 150, -400, 300] * u.arcsec,
                  [-500, 200, -250, 450, 100] * u.arcsec,
                  frame=aia_map.coordinate_frame)
five_point_crop = aia_map.crop(*points)

###############################################################################
# Plot both results with the coordinates used to define them.

fig = plt.figure(figsize=(10, 5))

ax1 = fig.add_subplot(121, projection=two_point_crop)
two_point_crop.plot(axes=ax1)
ax1.plot_coord(SkyCoord([bottom_left, top_right]), "o", color="cyan", fillstyle="none")
ax1.set_title("crop with 2 coordinates")

ax2 = fig.add_subplot(122, projection=five_point_crop)
five_point_crop.plot(axes=ax2)
ax2.plot_coord(points, "o", color="cyan", fillstyle="none")
ax2.set_title("crop with 5 coordinates")

plt.show()

###############################################################################
# In both cases, every input coordinate lies on or inside the result, which is
# the smallest pixel rectangle that contains them all.
