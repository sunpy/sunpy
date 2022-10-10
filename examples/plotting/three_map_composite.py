"""
=========================================
Creating a Composite Plot with Three Maps
=========================================

In this example, a composite plot is created with three maps.
It demonstrates how to specify contour levels, transparency, and
ordering when overlaying multiple maps.
"""
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
from sunpy.map import Map

###############################################################################
# First, we will import sample data from EIT, RHESSI, and AIA. The EIT data
# shows a hot region of the solar corona, while AIA shows the cooler upper
# region of the corona. RHESSI data is focused on a solar flare, and will be
# plotted using contours.
eit = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE)
rhessi = sunpy.map.Map(sunpy.data.sample.RHESSI_IMAGE)
aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

###############################################################################
# Next, specify the RHESSI contour levels to be plotted. When creating the
# plot, choose which Map will be used to create the WCS Axes that the other
# Maps will be plotted with. The EIT map is plotted first, followed by the
# RHESSI contours, and lastly the AIA map. Transparency is changed to 50% on
# the AIA map by specifying the parameter `alpha`, and the image data is
# autoaligned to the EIT WCS Axes. The parameter `zorder` specifies
# how each plot is layered (0 is plotted first and 1 is layered on top of 0,
# and so on). Lastly, note that `draw_grid()` is called using the EIT map.
levels = [50, 60, 70, 80, 90]*u.percent
fig = plt.figure()
ax = fig.add_subplot(projection=eit)
eit.plot(axes=ax, zorder=0)
rhessi.draw_contours(axes=ax, levels=levels, zorder=1)
aia.plot(axes=ax, alpha=0.5, autoalign=True, zorder=2)
eit.draw_grid(axes=ax)
plt.title("Three Map Composite")
plt.show()
