"""
====================================
Drawing a rotated rectangle on a map
====================================

This example will demonstrate how to draw a rectangle that can be rotated relative
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
center_coord = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=aia_map.coordinate_frame)
width = 400 * u.arcsec
height = 300 * u.arcsec

################################################################################
# To specify the rotation and center, we will use `~astropy.coordinates.SkyOffsetFrame`.
# Define the corner coordinates in the `~astropy.coordinates.SkyCoord` frame
# Transform corner coordinates to the original map's coordinate frame

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
    label='Rotated Rectangle'
)

plt.show()
