"""
===========================
Azimuthal Equidistant Projection with Reproject
===========================

In this example we use the `reproject` generate an image in `azimuthal equidistant projection <https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection>`__ from an AIA image.

You will need `reproject <https://reproject.readthedocs.io/en/stable/>`__ v0.6 or higher installed.
"""

import matplotlib.pyplot as plt
import reproject

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE

###############################################################################
# The SunPy Map factory accepts a wide variety of inputs for creating maps
# In this example we load the sample AIA image which we imported earlier from `sunpy.map.Map <https://docs.sunpy.org/en/stable/api/sunpy.map.Map.html>`
aia_map = sunpy.map.Map(AIA_171_IMAGE)

###############################################################################
# Negative value pixels can appear that lead to ugly looking images.
# This can be fix by setting the lower limit of normalization.
aia_map.plot_settings['norm'].vmin = 0
aia_map.plot_settings['norm'].vmax = 10000

###############################################################################
# Skycoord provides a flexible interface for celestial coordinate representation, manipulation, and transformation between systems.

# Here we define origin for the postel projection.
origin_hpc = SkyCoord(735*u.arcsec, -340*u.arcsec, frame=aia_map.coordinate_frame)
origin = origin_hpc.heliographic_stonyhurst

###############################################################################
# Create a FITS-WCS header from a coordinate object (SkyCoord) that is required to create a GenericMap.
out_shape = (750, 750)
out_header = sunpy.map.make_fitswcs_header(
    out_shape,
    origin,
    scale=[0.4, 0.4]*u.deg/u.pix,
    projection_code="ARC"
)
out_wcs = WCS(out_header)

###############################################################################
# Reproject the AIA image to the new WCS
output, _ = reproject.reproject_interp(aia_map, out_wcs, out_shape)

###############################################################################
# Create the reprojected Map
out_map = sunpy.map.Map(output, out_header)
out_map.plot_settings = aia_map.plot_settings

# Create a figure to plot the reprojected image
fig = plt.figure(figsize=(8, 4))


###############################################################################
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

# Show the plot
plt.show()
