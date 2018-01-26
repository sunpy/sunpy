"""
===============================
Fine grained Plotting Features
===============================

An example to show control over various plotting features.
"""
###############################################################################
# Import the necessary modules for plotting.

import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE


###############################################################################
# SkyCoord module provides flexible celestial coordinate representation and a
# draw_limb method overplots the boundary of the spacial objects. Date of the
# image taken can also be displayed in the plot.

aiamap = sunpy.map.Map(AIA_171_IMAGE)

bottom_left = SkyCoord(-400*u.arcsec, -900*u.arcsec, frame=aiamap.coordinate_frame)
top_right = SkyCoord(800*u.arcsec, 700*u.arcsec, frame=aiamap.coordinate_frame)
aiamap_sub = aiamap.submap(bottom_left, top_right)

title_obsdate = '{:%Y-%b-%d %H:%M:%S}'.format(aiamap_sub.date)

###############################################################################
# The figure displays coordinates on all four edges by default. The edges can be
# picked for labelling and in this example we remove the labels and tick on the
# top and right edges.
# The plots show the grid as white faint dashed lines by default and that can be removed.
# Also the coordinates are float values for general usage but integers can be used for either axis.
# For more information regarding the axis and grid settings, go to
# http://docs.astropy.org/en/stable/visualization/wcsaxes/ticks_labels_grid.html .

fig = plt.figure(figsize=(6, 6))
ax = plt.subplot(projection=aiamap_sub)
aiamap_sub.plot()
aiamap_sub.draw_limb(color='black', linewidth=2, linestyle='dashed')
overlay = ax.get_coords_overlay('heliographic_stonyhurst')
lon = overlay[0]
lat = overlay[1]

lon.set_ticks_visible(False)
lat.set_ticks_visible(False)
lat.set_ticklabel_visible(False)
lon.set_ticklabel_visible(False)

lon.coord_wrap = 180
lon.set_major_formatter('dd')

overlay.grid(color='blue', linewidth=2, linestyle='dashed')
ax.grid(False)
tx, ty = ax.coords
tx.set_major_formatter('s')
ty.set_major_formatter('s')
ax.set_title('AIA 171 $\AA$ {}'.format(title_obsdate))
ax.set_ylabel('Helioprojective Latitude [arcsec]')
ax.set_xlabel('Helioprojective Longitude [arcsec]')
plt.colorbar(fraction=0.045, pad=0.03, label='DN', ax=ax)
plt.show()
