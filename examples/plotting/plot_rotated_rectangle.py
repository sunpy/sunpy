"""
====================================
Drawing a rotated rectangle on a map
====================================

This example will demonstrate how to draw a rectangle that is rotated relative
to the axes on a map using :meth:`~sunpy.map.GenericMap.draw_quadrangle`.
"""
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame

import sunpy.data.sample
import sunpy.map

################################################################################
# Let's start with a sample image of AIA 171.

# sphinx_gallery_defer_figures

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

fig = plt.figure()
ax = fig.add_subplot(projection=aia_map)
aia_map.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)

################################################################################
# Define the rotation angle and center coordinate of the rectangle,
# as well as the width and height in physical units.

# sphinx_gallery_defer_figures

rotation_angle = 90 * u.deg
center_coord = SkyCoord(-300 * u.arcsec, -300 * u.arcsec, frame=aia_map.coordinate_frame)
width = 400 * u.arcsec
height = 300 * u.arcsec

################################################################################
# Now to specify the rotation and center for the rectangle, we will use `~astropy.coordinates.SkyOffsetFrame`.
# First we will define the bottom left and top right within the rotated frame
# and then transform these into the AIA 171 coordinate frame.

# sphinx_gallery_defer_figures

offset_frame = SkyOffsetFrame(origin=center_coord, rotation=rotation_angle)
corner_bottom_left = SkyCoord(lon=-width / 2, lat=-height / 2, frame=offset_frame)
corner_top_right = SkyCoord(lon=width / 2, lat=height / 2, frame=offset_frame)
corner_bottom_left = corner_bottom_left.transform_to(aia_map.coordinate_frame)
corner_top_right = corner_top_right.transform_to(aia_map.coordinate_frame)

################################################################################
# Finally, we will draw the rotated rectangle.

aia_map.draw_quadrangle(
    bottom_left=corner_bottom_left,
    top_right=corner_top_right,
    axes=ax,
    edgecolor="purple",
    linestyle="--",
    linewidth=2,
)

plt.show()
