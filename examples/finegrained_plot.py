"""
================
Fine grained Plotting Features
================

An example to show control over various plotting features.
"""
##############################################################################
# Import the necessary modules for plotting.

import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import figure
from astropy.coordinates import SkyCoord
import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE


###############################################################################
# skyCoord module provides flexible  celestial coordinate representation and a 
# draw_limb method overpplots the boundary of the spacial objects.Date of the 
# image taken can also be displayed in the plot. 

aiamap = sunpy.map.Map(AIA_171_IMAGE)

bl = SkyCoord(-400*u.arcsec, -900*u.arcsec, frame=aiamap.coordinate_frame)
tr = SkyCoord(800*u.arcsec, 700*u.arcsec, frame=aiamap.coordinate_frame)
ams = aiamap.submap(bl,tr)
title_obsdate='{:.20}'.format('{:%Y-%b-%d %H:%M:%s}'.format(ams.date))

###############################################################################
# The figure displays Coordinates on all four edges by default.The edges can be 
# picked for labelling and in this example we remove the labels and tick on the
# top and right edges.
# The plots show the grid as white faint dashed lines by default and that can be removed.
# Also the coordinates are float values for general usage but integers can be used for either axis.

fig = plt.figure(figsize=(5, 5))
ax = plt.subplot(projection=ams)
ams.plot()
ams.draw_limb(color='black',linewidth=2,linestyle='dashed')
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
ax.set_title('AIA 171$\AA$ '+ title_obsdate)
ax.set_ylabel('y [arcsec]')
ax.set_xlabel('x [arcsec]')
plt.colorbar(fraction=0.035, pad=0.03,label='DN', ax=ax)
# plt.colorbar(ax=ax)
plt.show()



