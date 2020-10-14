"""
==========================================================
Drawing a Rectangle on a map
==========================================================

In this example we are going to look at how we can plot rectangles on maps
while describing the dimensions using different methods.
"""
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map

################################################################################
# We use the `~sunpy.map.GenericMap.draw_rectangle` method to draw the rectangles.
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
fig = plt.figure(figsize=(5, 5))
aia_map.plot()

################################################################################
# The first option is to provide the top left and bottom right corners as a
# single `~astropy.coordinates.SkyCoord`.
fig.add_subplot(111, projection=aia_map)
coords = SkyCoord(
    Tx=(-400, 300) * u.arcsec,
    Ty=(-300, 200) * u.arcsec,
    frame=aia_map.coordinate_frame,
)
aia_map.draw_rectangle(
    coords,
    color="blue",
    linestyle="-",
)
aia_map.plot()

################################################################################
# The second option is to pass the bottom left and top right corners as
# separate `~astropy.coordinates.SkyCoord` objects.
fig.add_subplot(111, projection=aia_map)
bottom_left = SkyCoord(-400 * u.arcsec, -300 * u.arcsec, frame=aia_map.coordinate_frame)
top_right = SkyCoord(300 * u.arcsec, 200 * u.arcsec, frame=aia_map.coordinate_frame)
aia_map.draw_rectangle(
    bottom_left,
    top_right=top_right,
    color="green",
    linestyle="--",
)
aia_map.plot()

################################################################################
# The third option is to provide the coordinate of the bottom left corner along with the
# width and the height of the rectangle.
fig.add_subplot(111, projection=aia_map)
bottom_left = SkyCoord(-400 * u.arcsec, -300 * u.arcsec, frame=aia_map.coordinate_frame)
width = 700 * u.arcsec
height = 500 * u.arcsec
aia_map.draw_rectangle(
    bottom_left,
    width=width,
    height=height,
    color="yellow",
    linestyle="-.",
)
aia_map.plot()

################################################################################
# We also have a fourth option - to specify the dimensions in terms of pixel coordinates.
#
# To do this, we use the `~sunpy.map.GenericMap.pixel_to_world` method to convert the
# pixel coordinates to world coordinates.
fig.add_subplot(111, projection=aia_map)
bottom_left = aia_map.pixel_to_world(350 * u.pixel, 390 * u.pixel)
top_right = aia_map.pixel_to_world(650 * u.pixel, 600 * u.pixel)
aia_map.draw_rectangle(
    bottom_left,
    top_right=top_right,
    color="red",
    linestyle=":",
)
aia_map.plot()
