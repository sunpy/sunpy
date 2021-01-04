"""
==================================
Set Axis Range When Plotting a Map
==================================

In this example we are going to look at how to set the axes
range using Matplotlib's `set_xlim` and `set_ylim` when plotting a
Map with WCSAxes.
"""

import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map

###############################################################################
# Lets start by creating a Map from the sample data.
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

###############################################################################
# Now lets say for example we are only interested in plotting a certain region
# of this Map. One way this could be done is to create a submap over the region
# of interest and then plotting that. Another useful way is to set the axes
# range over which to plot using Matplotlib's `Axes.set_xlim` and `Axes.set_ylim` functionality.
# The axes that Matplotlib uses is in pixel coordinates (e.g. of image data array)
# rather than world coordinates (e.g. in arcsecs) so we need to define our limits that
# are passed to set_xlim(), set_lim() to pixel coordinates.
# We can define our limits we want to use in world coordinates and then work out what pixel
# coordinates these correspond to.
# Lets choose xlimits and ylimits in arcsecs that we are interested in.
xlims_world = [500, 1100]*u.arcsec
ylims_world = [-800, 0]*u.arcsec

###############################################################################
# We can then covert these into a SkyCoord which can be passed to :func:`~sunpy.map.GenericMap.world_to_pixel` to
# determine which pixel coordinates these represent on the Map.
world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_map.coordinate_frame)
pixel_coords = aia_map.world_to_pixel(world_coords)

# we can then pull out the x and y values of these limits.
xlims_pixel = pixel_coords.x.value
ylims_pixel = pixel_coords.y.value

###############################################################################
# We can now plot this Map and then use the x_lims_pixel and y_lims_pixel to set
# the range of the axes for which to plot.
fig = plt.figure()
ax = plt.subplot(projection=aia_map)
aia_map.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

ax.set_xlim(xlims_pixel)
ax.set_ylim(ylims_pixel)
plt.show()
