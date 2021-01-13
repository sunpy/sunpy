"""
=====================================
Plotting points on a Map with WCSAxes
=====================================

This example demonstrates the plotting of a point, a line and an arch in pixel coordinates,
world coordinates, and SkyCoords respectively when plotting a map with WCSAxes.
"""
import matplotlib.pyplot as plt
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map
from sunpy.coordinates.utils import GreatArc

###############################################################################
# We will start by creating a map using an AIA 171 image.
# Now we will plot a line on the map by using coordinates in arcseconds.
# The array below `xx` and `yy` are the x and y coordinates that define a
# line from the Sun center (at 0, 0) to the point (500, 500) in arcsecs.
# Plotting with WCSAxes expects pixel coordinates, but we can plot
# world coordinates (i.e. arcsec) by using the ``transform`` keyword.
my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection=my_map)
my_map.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

xx = np.arange(0, 500)
yy = xx

# Note that the coordinates need to be in degrees rather than arcseconds.
ax.plot(xx*u.arcsec.to(u.deg), yy*u.arcsec.to(u.deg),
        color='r',
        transform=ax.get_transform("world"),
        label=f'WCS coordinate [{0*u.arcsec}, {500*u.arcsec}]')

# Here we will plot a point in pixel coordinates (i.e. array index).
# We are defining a pixel coordinate in the middle of the image.
pixel_coord = [my_map.data.shape[0]/2., my_map.data.shape[1]/2.] * u.pix
ax.plot(pixel_coord[0], pixel_coord[1], 'x', color='w',
        label=f'Pixel coordinate [{pixel_coord[0]}, {pixel_coord[1]}]')
ax.coords.grid(color='yellow', linestyle='solid', alpha=0.5)

# Now let's plot a point and an arc on map using two separate SkyCoords.
# This will plot a point (at -250,-250) on the map using a SkyCoord.
ax.plot_coord(SkyCoord(-250*u.arcsec, -250*u.arcsec, frame=my_map.coordinate_frame), "o",
              label=f'SkyCoord [{-250*u.arcsec}, {-250*u.arcsec}]')

# Finally, this will plot an arc using a SkyCoord.
start = SkyCoord(723 * u.arcsec, -500 * u.arcsec, frame=my_map.coordinate_frame)
end = SkyCoord(-100 * u.arcsec, 900 * u.arcsec, frame=my_map.coordinate_frame)

great_arc = GreatArc(start, end)

my_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
ax.plot_coord(great_arc.coordinates(), color='c',
              label=f'SkyCoord [{723*u.arcsec}, {-500*u.arcsec}],\n \
                               [{-100*u.arcsec}, {900*u.arcsec}]')

plt.legend(loc="lower center")
plt.show()
