"""
========================
Creating Carrington Maps
========================

In this example we use the `reproject` generate a map in heliographic Carrington coordinates from a full-disk AIA image.

You will need `reproject <https://reproject.readthedocs.io/en/stable/>`__ v0.6 or higher installed.
"""
# sphinx_gallery_thumbnail_number = 2

import matplotlib.pyplot as plt

import sunpy.data.sample
import sunpy.map
from sunpy.map.header_helper import make_heliographic_header

###############################################################################
# We will start with using sunpy's sample data for this example.

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE)

fig = plt.figure()
ax = fig.add_subplot(projection=aia_map)
aia_map.plot(axes=ax)

###############################################################################
# Reproject works by transforming an input image to a desired World Coordinate
# System (WCS) projection. Here we use :func:`sunpy.map.header_helper.make_heliographic_header`
# to create a FITS WCS header based on a heliographic Carrington reference
# coordinate.

shape = (720, 1440)
carr_header = make_heliographic_header(aia_map.date, aia_map.observer_coordinate, shape, frame='carrington')

###############################################################################
# With the new header, re-project the data into the new coordinate system.
# The :meth:`~sunpy.map.GenericMap.reproject_to` defaults to using
# the fast :func:`reproject.reproject_interp` algorithm, but a different
# algorithm can be specified (e.g., :func:`reproject.reproject_adaptive`).

outmap = aia_map.reproject_to(carr_header)

###############################################################################
# Plot the result.

fig = plt.figure()
ax = fig.add_subplot(projection=outmap)
outmap.plot(axes=ax)
outmap.draw_limb(color='blue')

plt.show()
