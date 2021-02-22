"""
============================
Drawing a rectangle on a map
============================

This example will demonstrate how to draw a rectangle on a map using :meth:`~sunpy.map.GenericMap.draw_rectangle`.
"""
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map

################################################################################
# Let's start with a sample AIA image.

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

################################################################################
# Here are four different ways to draw a rectangle. The first three ways
# directly calls the `~astropy.coordinates.SkyCoord` class. The fourth way
# converts pixel coordinates to the equivalent `~astropy.coordinates.SkyCoord`
# objects using the :meth:`~sunpy.map.GenericMap.pixel_to_world`.

fig = plt.figure(figsize=(5, 5))
fig.add_subplot(111, projection=aia_map)
aia_map.plot(clip_interval=(1, 99.99)*u.percent)

# Specify two opposite corners of the rectangle as a single, two-element
# SkyCoord object
coords = SkyCoord(
    Tx=(100, 500) * u.arcsec,
    Ty=(200, 500) * u.arcsec,
    frame=aia_map.coordinate_frame,
)
aia_map.draw_rectangle(
    coords,
    edgecolor="blue",
    linestyle="-",
    linewidth=2,
    label='2-element SkyCoord'
)

# Specify two opposite corners of the rectangle as separate SkyCoord objects
bottom_left = SkyCoord(-500 * u.arcsec, 200 * u.arcsec, frame=aia_map.coordinate_frame)
top_right = SkyCoord(-100 * u.arcsec, 500 * u.arcsec, frame=aia_map.coordinate_frame)
aia_map.draw_rectangle(
    bottom_left,
    top_right=top_right,
    edgecolor="green",
    linestyle="--",
    linewidth=2,
    label='two SkyCoords'
)

# Specify one corner of the rectangle and the rectangle's width and height
bottom_left = SkyCoord(-500 * u.arcsec, -500 * u.arcsec, frame=aia_map.coordinate_frame)
width = 400 * u.arcsec
height = 300 * u.arcsec
aia_map.draw_rectangle(
    bottom_left,
    width=width,
    height=height,
    edgecolor="yellow",
    linestyle="-.",
    linewidth=2,
    label='width/height'
)

# Draw a desired rectangle in pixel coordinates by first converting to SkyCoord objects
bottom_left = aia_map.pixel_to_world(600 * u.pixel, 350 * u.pixel)
top_right = aia_map.pixel_to_world(800 * u.pixel, 450 * u.pixel)
aia_map.draw_rectangle(
    bottom_left,
    top_right=top_right,
    edgecolor="red",
    linestyle=":",
    linewidth=2,
    label='pixel_to_world()'
)

plt.legend()
plt.show()
