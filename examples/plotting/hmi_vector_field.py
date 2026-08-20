"""
==============================================================
Overlaying the HMI vector field on a line-of-sight magnetogram
==============================================================

In this example we will download the vector field data from a Spaceweather HMI
Active Region Patch (SHARP) dataset and overlay it as vectors
on top of the line-of-sight magnetogram.
"""
# sphinx_gallery_tags = ["Map", "HMI", "SHARP", "JSOC", "Visualization", "Active Regions"]

import os

import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

##############################################################################
# To get the vector field data, we need to query the
# `JSOC <http://jsoc.stanford.edu/>`__ using `Fido
# <sunpy.net.fido_factory.UnifiedDownloaderFactory>` and the search
# attributes in `sunpy.net.jsoc`.
#
# The ``hmi.sharp_720s`` series provides the full magnetic field vector as
# three segments: the field strength ("field"), inclination and azimuth,
# alongside the line-of-sight magnetogram ("magnetogram").
#
# Here we download all four for HARP number 4536.
#
# Exporting data from the JSOC requires registering your email first.
# Please replace this with your email address once you have registered
# like so: jsoc_email = "your_email@example.com"
# See `this page <http://jsoc.stanford.edu/ajax/register_email.html>`__ for more details.

jsoc_email = os.environ["JSOC_EMAIL"]

result = Fido.search(
    a.Time("2014-09-10 17:00:00", "2014-09-10 17:01:00"),
    a.jsoc.Series("hmi.sharp_720s"),
    a.jsoc.PrimeKey("HARPNUM", 4536),
    a.jsoc.Notify(jsoc_email),
    a.jsoc.Segment("magnetogram") & a.jsoc.Segment("field")
    & a.jsoc.Segment("inclination") & a.jsoc.Segment("azimuth"),
)
print(result)

##############################################################################
# Now we can download the files and create a few maps.

files = Fido.fetch(result)
# Sorting the filenames puts the segments in alphabetical order:
# azimuth, field, inclination, magnetogram.
azimuth_map, field_map, inclination_map, magnetogram_map = sunpy.map.Map(sorted(files))

##############################################################################
# The magnetic field vector is described in a spherical form: the field
# strength, its inclination to the line of sight, and the azimuth of its
# transverse component, measured counterclockwise from the "up" direction
# (+y) of the CCD image. The SHARP azimuth already has its 180-degree
# ambiguity resolved :cite:p:`hoeksema_helioseismic_2014`. Equation 1 of
# :cite:t:`sun_coordinate_2013` gives the transverse components along the
# CCD axes:

inclination = np.deg2rad(inclination_map.data)
azimuth = np.deg2rad(azimuth_map.data)

b_x = -field_map.data * np.sin(inclination) * np.sin(azimuth)
b_y = field_map.data * np.sin(inclination) * np.cos(azimuth)

##############################################################################
# HMI images have solar north pointing approximately down, so we rotate the
# maps to put north up. We rotate the component arrays rather than the
# azimuth map, as interpolating angles across the 0/360 degree wrap would
# corrupt them; reusing the magnetogram's metadata supplies the coordinate
# information. Rotating a map only resamples its values, so we must also
# rotate the vectors into the new image axes using the map's
# :attr:`~sunpy.map.GenericMap.rotation_matrix`, which for HMI is nearly
# (but not exactly) a 180-degree rotation, i.e. a sign flip of both
# components.

rmatrix = magnetogram_map.rotation_matrix
b_x_ccd = sunpy.map.Map(b_x, magnetogram_map.meta).rotate().data
b_y_ccd = sunpy.map.Map(b_y, magnetogram_map.meta).rotate().data
b_x = rmatrix[0, 0] * b_x_ccd + rmatrix[0, 1] * b_y_ccd
b_y = rmatrix[1, 0] * b_x_ccd + rmatrix[1, 1] * b_y_ccd
magnetogram_map = magnetogram_map.rotate()

##############################################################################
# Plotting an arrow at every pixel would be unreadable, so we sample the
# field on a coarser grid. We also mask out weak transverse fields, which
# are dominated by noise, and the strongest fields, so that a few long
# arrows do not dominate the plot.

ny, nx = magnetogram_map.data.shape
step = 5

yy, xx = np.mgrid[0:ny:step, 0:nx:step]
b_x = b_x[::step, ::step]
b_y = b_y[::step, ::step]
b_transverse = np.hypot(b_x, b_y)

b_min = 300
b_max = 1500
good = (b_transverse > b_min) & (b_transverse < b_max)

##############################################################################
# Finally, we plot the line-of-sight magnetogram and overlay the transverse
# field as arrows. On a map projection, ``quiver`` works in pixel
# coordinates by default, which is what our sampled grid is in.

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection=magnetogram_map)
# Adjust the vmin/vmax of the normalization like so:
magnetogram_map.plot_settings["norm"].vmin = -1500
magnetogram_map.plot_settings["norm"].vmax = 1500
magnetogram_map.plot(axes=ax, cmap="Greys_r")
ax.set_title(r'$B_{LOS}$ & Vector Field - ' + magnetogram_map.date.isot)

# We want to crop the plot to a tighter region of interest
corners = SkyCoord(Tx=[-150, 0] * u.arcsec, Ty=[80, 230] * u.arcsec,
                   frame=magnetogram_map.coordinate_frame)
x_pix, y_pix = magnetogram_map.wcs.world_to_pixel(corners)
ax.set_xlim(x_pix)
ax.set_ylim(y_pix)

ax.quiver(xx[good], yy[good], b_x[good], b_y[good], color="red", alpha=0.5, pivot="middle")

plt.show()
