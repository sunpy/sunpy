"""
====================
HMI Showcase: Cutout
====================

This example demonstrates how to plot a cutout region of a `~sunpy.map.Map` with connector lines
that indicate the region of interest in the full-disk image.

Since this example deals with the creation of a specific style of image, there are multiple lines that deal directly with matplotlib axes.
Furthermore, unlike other examples, since this deals with the creation of one figure, there are fewer breaks in the example.
"""
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.data.sample import HMI_LOS_IMAGE

##############################################################################
# Firstly, we will use the HMI LOS image within the sunpy sample data and focus
# the cutout over an active region near the solar center.

magnetogram = sunpy.map.Map(HMI_LOS_IMAGE).rotate()
left_corner = SkyCoord(Tx=-142*u.arcsec, Ty=50*u.arcsec, frame=magnetogram.coordinate_frame)
right_corner = SkyCoord(Tx=158*u.arcsec, Ty=350*u.arcsec, frame=magnetogram.coordinate_frame)

##############################################################################
# As there is "data" off the limb of the magetogram, we will need to mask it away.
# We will create a second sunpy map that holds this masked data for use later.
# From here the rest of the comments will be inside the code block.

x, y = np.meshgrid(*[np.arange(v.value) for v in magnetogram.dimensions]) * u.pixel
hpc_coords = magnetogram.pixel_to_world(x, y)
solar_rad = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / magnetogram.rsun_obs
mask = np.ma.masked_greater(solar_rad, 1)
magnetogram_big = sunpy.map.Map(magnetogram.data, magnetogram.meta, mask=mask.mask)


# The next stage is to setup the figure in two stages.
# The first stage is dealing with the full disk image.

fig = plt.figure()

# We want to create a nice normalation range for the image

norm = matplotlib.colors.SymLogNorm(50, vmin=-7.5e2, vmax=7.5e2)

# Plot the first magnetogram

ax1 = fig.add_subplot(121, projection=magnetogram_big)
magnetogram_big.plot(axes=ax1, cmap='RdBu_r', norm=norm, annotate=False,)

# These lines deal with hiding the axis, its ticks and labels

lon, lat = ax1.coords[0], ax1.coords[1]
lon.frame.set_linewidth(0)
lat.frame.set_linewidth(0)
lon.set_ticks_visible(False)
lat.set_ticks_visible(False)
lon.set_ticklabel_visible(False)
lat.set_ticklabel_visible(False)

# We will draw the rectangle around the region we plan to showcase in the cutout image.

magnetogram_big.draw_rectangle(left_corner, height=right_corner.Tx-left_corner.Tx,
                               width=right_corner.Ty-left_corner.Ty, color='k', lw=1)
# As well as draw the full cooridate grid on top.

magnetogram_big.draw_grid(axes=ax1, color='k', alpha=0.25, lw=0.5)


# Following this, we will now deal with the zoomed-in magnetogram.

magnetogram_small = magnetogram.submap(left_corner, top_right=right_corner)
ax2 = fig.add_subplot(122, projection=magnetogram_small)
im = magnetogram_small.plot(axes=ax2, norm=norm, cmap='RdBu_r', annotate=False,)
ax2.grid(alpha=0)
lon, lat = ax2.coords[0], ax2.coords[1]

# Unlike the full disk image, here we will just clean up the axis labels and ticks.

lon.frame.set_linewidth(1)
lat.frame.set_linewidth(1)
lon.set_ticklabel()
lat.set_ticklabel(rotation='vertical',)
lon.set_axislabel('Helioprojective Longitude',)
lat.set_axislabel('Helioprojective Latitude',)
lat.set_axislabel_position('r')
lat.set_ticks_position('r')
lat.set_ticklabel_position('r')

# Now for the finishing touches, we want to add two lines that will connect
# the two images as well as a colorbar.

xpix, ypix = magnetogram_big.world_to_pixel(right_corner)
con1 = ConnectionPatch(
    (0, 1), (xpix.value, ypix.value), 'axes fraction', 'data', axesA=ax2, axesB=ax1,
    arrowstyle='-', color='k', lw=1
)
xpix, ypix = magnetogram_big.world_to_pixel(
    SkyCoord(right_corner.Tx, left_corner.Ty, frame=magnetogram_big.coordinate_frame))
con2 = ConnectionPatch(
    (0, 0), (xpix.value, ypix.value), 'axes fraction', 'data', axesA=ax2, axesB=ax1,
    arrowstyle='-', color='k', lw=1
)
ax2.add_artist(con1)
ax2.add_artist(con2)

pos = ax2.get_position().get_points()
cax = fig.add_axes([
    pos[0, 0], pos[1, 1]+0.01, pos[1, 0]-pos[0, 0], 0.025
])
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', label="abc")

# For the colorbar we want it to have three fixed ticks

cbar.locator = matplotlib.ticker.FixedLocator([-1e2, 0, 1e2])
cbar.set_label("LOS Magnetic Field [Gauss]", labelpad=-40, rotation=0)
cbar.update_ticks()
cbar.ax.xaxis.set_ticks_position('top')

plt.show()
