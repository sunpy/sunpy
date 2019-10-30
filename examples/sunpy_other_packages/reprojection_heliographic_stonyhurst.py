"""
===========================
Creating a Heliographic Map
===========================

In this example we use the `reproject` package to make a heliographic image
from an AIA image.

You will need to install reproject to run this example.
"""
# sphinx_gallery_thumbnail_number = 2

import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from reproject import reproject_interp

import sunpy.data.sample
import sunpy.map

###############################################################################
# We start with the sample data
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE)

fig = plt.figure()
ax = plt.subplot(projection=aia_map)
aia_map.plot(ax)

###############################################################################
# Reproject works by transforming an input image (with a `~astropy.wcs.WCS`) to
# a output image, specified by a different WCS object. Therefore we need to
# build a `~astropy.wcs.WCS` object describing the output we desire.
# To do this we use the `sunpy.map.make_fitswcs_header` which assists us in
# constructing this WCS object. Here we create a WCS based on a heliographic
# Stonyhurst reference coordinate and with the CAR (plate carrée) projection.
shape_out = [720, 1440]
header = sunpy.map.make_fitswcs_header(shape_out,
                                       SkyCoord(0, 0, unit=u.deg,
                                                frame="heliographic_stonyhurst",
                                                obstime=aia_map.date),
                                       scale=[180 / shape_out[0],
                                              360 / shape_out[1]] * u.deg / u.pix,
                                       projection_code="CAR")

out_wcs = WCS(header)

###############################################################################
# With the new header, re-project the data into the new coordinate system.
# Here we are using the fastest but least accurate method of reprojection,
# `reproject.reproject_interp`, a more accurate but slower method is
# `reproject.reproject_adaptive`.
array, footprint = reproject_interp(aia_map, out_wcs, shape_out=shape_out)
outmap = sunpy.map.Map((array, header))
outmap.plot_settings = aia_map.plot_settings

###############################################################################
# Plot the result
fig = plt.figure()
ax = plt.subplot(projection=outmap)
outmap.plot(ax)

ax.set_xlim(0, shape_out[1])
ax.set_ylim(0, shape_out[0])

plt.show()
