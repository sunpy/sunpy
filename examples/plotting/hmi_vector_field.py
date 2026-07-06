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
# Now we can download the files and load them into maps.

files = Fido.fetch(result)
# Sorting the filenames puts the segments in alphabetical order:
# azimuth, field, inclination, magnetogram.
azimuth_map, field_map, inclination_map, magnetogram_map = sunpy.map.Map(sorted(files))

##############################################################################
# The field strength, inclination and azimuth describe the magnetic field
# vector in spherical coordinates, so we first convert it to a Cartesian
# form. The azimuth is measured counterclockwise from the "up" direction
# (+y) of the CCD image, so the components of the field transverse to the
# line of sight are given by (see equation 1 of :cite:t:`sun_coordinate_2013`):

inclination = np.deg2rad(inclination_map.data)
azimuth = np.deg2rad(azimuth_map.data)

b_x = -field_map.data * np.sin(inclination) * np.sin(azimuth)
b_y = field_map.data * np.sin(inclination) * np.cos(azimuth)

##############################################################################
# HMI images have solar north pointing down, so we rotate the maps to put
# solar north up. We deliberately computed the vector components *before* rotating
# them, so that the sign of both components is flipped when we rotate the
# magnetogram map. The azimuth map directly would interpolate its values
# across the 0/360 degree wrap and corrupt them. Instead, we place the components
# into maps that share the magnetogram's coordinate information and rotate those.
# Rotating the image by 180 degrees also rotates the direction each vector points,
# which amounts to flipping the sign of both components.

b_x = -sunpy.map.Map(b_x, magnetogram_map.meta).rotate().data
b_y = -sunpy.map.Map(b_y, magnetogram_map.meta).rotate().data
magnetogram_map = magnetogram_map.rotate()

##############################################################################
# Plotting an arrow at every pixel would be unreadable, so we sample the
# field on a coarser grid. We also mask out weak transverse fields (which
# are dominated by noise) as well as unphysically strong ones.

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
# Adjust the vmin/vmax of the normalisation like so:
magnetogram_map.plot_settings["norm"].vmin = -1500
magnetogram_map.plot_settings["norm"].vmax = 1500
magnetogram_map.plot(axes=ax, cmap="Greys_r")
ax.set_title(r'$B_{LOS}$ & Vector Field - ' + magnetogram_map.date.isot)

# We want to crop the plot to a tighter region of interest.
corners = SkyCoord(Tx=[-150, 0] * u.arcsec, Ty=[80, 230] * u.arcsec,
                   frame=magnetogram_map.coordinate_frame)
x_pix, y_pix = magnetogram_map.wcs.world_to_pixel(corners)
ax.set_xlim(x_pix)
ax.set_ylim(y_pix)

ax.quiver(xx[good], yy[good], b_x[good], b_y[good], color="red", alpha=0.5)

plt.show()
