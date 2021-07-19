"""
===========================
Creating a Heliographic Map
===========================

In this example we use the `reproject` generate an image in heliographic coordinates from an AIA image.

You will need `reproject <https://reproject.readthedocs.io/en/stable/>`__ v0.6 or higher installed.
"""
# sphinx_gallery_thumbnail_number = 2

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map

###############################################################################
# We will start with using sunpy's sample data for this example.

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE)

fig = plt.figure()
ax = plt.subplot(projection=aia_map)
aia_map.plot(ax)

###############################################################################
# Reproject works by transforming an input image to a desired World Coordinate
# System (WCS) projection. Here we use :func:`sunpy.map.make_fitswcs_header`
# to create a FITS WCS header based on a Stonyhurst heliographic reference
# coordinate and the CAR (plate carrée) projection.

shape_out = (720, 1440)
frame_out = SkyCoord(0, 0, unit=u.deg,
                     frame="heliographic_stonyhurst",
                     obstime=aia_map.date,
                     rsun=aia_map.coordinate_frame.rsun)
header = sunpy.map.make_fitswcs_header(shape_out,
                                       frame_out,
                                       scale=(360 / shape_out[1],
                                              180 / shape_out[0]) * u.deg / u.pix,
                                       projection_code="CAR")

###############################################################################
# With the new header, re-project the data into the new coordinate system.
# The :meth:`~sunpy.map.GenericMap.reproject_to` defaults to using
# the fast :func:`reproject.reproject_interp` algorithm, but a different
# algorithm can be specified (e.g., :func:`reproject.reproject_adaptive`).

outmap = aia_map.reproject_to(header)

###############################################################################
# Plot the result.

fig = plt.figure()
ax = plt.subplot(projection=outmap)
outmap.plot(ax)
outmap.draw_limb(color='blue')

plt.show()
