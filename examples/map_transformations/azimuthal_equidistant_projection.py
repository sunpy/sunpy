"""
===============================================
Azimuthal Equidistant Projection with Reproject
===============================================

In this example, we will reproject an AIA image into an
`azimuthal equidistant projection <https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection>`__.
You will need v0.6 or higher  of the `reproject` package installed.

"""

import matplotlib.pyplot as plt
import reproject

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE

###############################################################################
# `~sunpy.map.Map` accepts a wide variety of inputs.
# In this example we load the sample AIA image which we imported earlier from
# `sunpy.data.sample`.

aia_map = sunpy.map.Map(AIA_171_IMAGE)

###############################################################################
# Negative value pixels can appear that lead to ugly looking images.
# This can be fixed by setting the lower limit of the normalization.

aia_map.plot_settings['norm'].vmin = 0
aia_map.plot_settings['norm'].vmax = 10000

###############################################################################
# Next, we define the origin of the postel projection using an `astropy.coordinates.SkyCoord` object.
# `~astropy.coordinates.Skycoord` provides a flexible interface for celestial coordinate representation, manipulation, and transformation between systems.

origin_hpc = SkyCoord(735*u.arcsec, -340*u.arcsec, frame=aia_map.coordinate_frame)
origin = origin_hpc.heliographic_stonyhurst

###############################################################################
# We can then use our origin coordinate to create a FITS-WCS header.
# From this new header, we can derive the world coordinate system (WCS)
# of our postel projection.

out_shape = (750, 750)
out_header = sunpy.map.make_fitswcs_header(
    out_shape,
    origin,
    scale=[0.4, 0.4]*u.deg/u.pix,
    projection_code="ARC"
)
out_wcs = WCS(out_header)

###############################################################################
# Now, we can use the WCS of postel projection to reproject the AIA map.

out_map = aia_map.reproject_to(out_wcs)
out_map.plot_settings = aia_map.plot_settings

###############################################################################

# sphinx_gallery_defer_figures

# Finally, we'll plot both our original and reprojected images.

fig = plt.figure(figsize=(8, 4))

###############################################################################

# sphinx_gallery_defer_figures

# Plot the original AIA map

ax = fig.add_subplot(1, 2, 1, projection=aia_map)
aia_map.plot(axes=ax)
aia_map.draw_grid(axes=ax, color='blue')
aia_map.draw_limb(axes=ax, color='blue')
ax.plot_coord(origin, 'o', color='red', fillstyle='none', markersize=20)

###############################################################################
# Plot the reprojected AIA map

ax = fig.add_subplot(1, 2, 2, projection=out_map)
out_map.plot(axes=ax, autoalign=True, annotate=False)
out_map.draw_grid(axes=ax, color='blue')
out_map.draw_limb(axes=ax, color='blue')
ax.set_xlim(-0.5, out_shape[1] - 0.5)
ax.set_ylim(-0.5, out_shape[0] - 0.5)
ax.set_title('Postel projection centered at ROI', y=-0.1)
plt.show()
