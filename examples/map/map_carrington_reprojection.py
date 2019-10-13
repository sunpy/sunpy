"""
==================
Re-projecting Maps
==================

How to reproject a map into
"""
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from reproject import reproject_interp

import sunpy.map
import sunpy.data.sample

###############################################################################
# We start with the sample data
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE)

fig = plt.figure()
ax = plt.subplot(projection=aia_map)
aia_map.plot(ax)

###############################################################################
# Now construct a header for the output
shape_out = [720, 1440]
header = sunpy.map.make_fitswcs_header(np.empty(shape_out),
                                       SkyCoord(0, 0, unit=u.deg,
                                                frame="heliographic_stonyhurst",
                                                obstime=aia_map.date),
                                       scale=[180 / shape_out[0],
                                              360 / shape_out[1]] * u.deg / u.pix,
                                       projection_code="CAR")

out_wcs = WCS(header)

###############################################################################
# With the new header, re-project the data into the new coordinate system.
array, footprint = reproject_interp(
    (aia_map.data, aia_map.wcs), out_wcs, shape_out=shape_out)
outmap = sunpy.map.Map((array, header))
outmap.plot_settings = aia_map.plot_settings

###############################################################################
# Plot the result
fig = plt.figure()
ax = plt.subplot(projection=outmap.wcs)
outmap.plot(ax)

ax.set_xlim(0, shape_out[1])
ax.set_ylim(0, shape_out[0])

plt.show()
