"""
=========================
Zoomed-in Cutout of a Map
=========================

This example demonstrates how to plot a cutout region of a `~sunpy.map.Map` with
connector lines that indicate the region of interest in the full-disk image.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE

###############################################################################
# We will start off by defining a region of interest that we want to focus on.
# In this case, we want to focus on the Coronal Mass Ejection (CME) in the lower right.
length = 600 * u.arcsec
x0 = 550 * u.arcsec
y0 = -700 * u.arcsec

# Now we will create a Map and then crop the region of interest.
smap = sunpy.map.Map(AIA_171_IMAGE)
bottom_left = SkyCoord(x0, y0, frame=smap.coordinate_frame)
top_right = SkyCoord(x0 + length, y0+length,
                     frame=smap.coordinate_frame)

submap = smap.submap(bottom_left, top_right=top_right)

###############################################################################
# Now we have the submap of the region of interest. We can focus on creating
# the full image.
fig = plt.figure(figsize=(18, 7))
ax1 = fig.add_subplot(1, 2, 1, projection=smap)

smap.plot(axes=ax1, clip_interval=(1, 99.9)*u.percent)

# Draw a box on the image highlighting the region of interest.
smap.draw_rectangle(bottom_left, color='k', height=length, width=length)

# Create a second axis on the plot which will be used to create the zoomed-in image.
ax2 = fig.add_subplot(1, 2, 2, projection=submap)

submap.plot(axes=ax2, clip_interval=(1, 99.9)*u.percent)
submap.draw_grid(grid_spacing=20*u.deg)

# We want to remove some duplicate text on the cropped image.
ax2.set_title(None)
ax2.set_ylabel('Solar Latitude')
ax2.set_xlabel('', color='w')

# The final touch is the creation of the connecting lines between the two images.
xpix, ypix = smap.world_to_pixel(top_right)
con1 = ConnectionPatch(
    (0, 1), (xpix.value, ypix.value), 'axes fraction', 'data', axesA=ax2, axesB=ax1,
    arrowstyle='-', color='k', lw=1
)
xpix, ypix = smap.world_to_pixel(SkyCoord(top_right.Tx, bottom_left.Ty,
                                          frame=smap.coordinate_frame))
con2 = ConnectionPatch(
    (0, 0), (xpix.value, ypix.value), 'axes fraction', 'data', axesA=ax2, axesB=ax1,
    arrowstyle='-', color='k', lw=1
)
ax2.add_artist(con1)
ax2.add_artist(con2)
plt.show()
