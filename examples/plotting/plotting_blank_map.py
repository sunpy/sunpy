"""
=============================
Plot positions on a blank map
=============================

This example showcases how to plot positions on a blank map.
It is often useful to plot coordinate positions of events on a blank helioprojective coordinate map.
In this example, we create an empty map with a WCS defined by a helioprojective frame as observed from
Earth at a certain time, and show how you can plot different coordinates on it.
"""
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.coordinates import frames

################################################################################
# First we will create a blank map using with an array of zeros.
# Since there is no WCS information, we will need to construct a header to pass to Map.

data = np.full((10, 10), np.nan)

# Define a reference coordinate and create a header using sunpy.map.make_fitswcs_header
skycoord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime='2013-10-28',
                    observer='earth', frame=frames.Helioprojective)

# Scale set to the following for solar limb to be in the field of view
header = sunpy.map.make_fitswcs_header(data, skycoord, scale=[220, 220]*u.arcsec/u.pixel)

# Use sunpy.map.Map to create the blank map
blank_map = sunpy.map.Map(data, header)

################################################################################
# Now we have constructed the map, we can plot it and mark important locations to it.
# Initialize the plot and add the map to it

# sphinx_gallery_defer_figures

fig = plt.figure()
ax = fig.add_subplot(projection=blank_map)
blank_map.plot(axes=ax)
blank_map.draw_limb(axes=ax, color="k")
blank_map.draw_grid(axes=ax, color="k")

################################################################################
# Coordinates that are being plotted - (0, 0), (50, 100) and (400, 400).

# sphinx_gallery_defer_figures

xc = [0, 50, 400] * u.arcsec
yc = [0, 100, 400] * u.arcsec

################################################################################
# Note that the blue dot corresponding to (0, 0) in helioprojective coordinates is
# not at the intersection of the heliographic equator and the heliographic prime
# meridian. (0, 0) in helioprojective coordinates *is* at the center of the solar
# disk, and at the center of the overall image. The reason that the intersection of
# the two heliographic grid lines is not also at disk center is because the observer
# (specified as Earth in this example) is at non-zero heliographic latitude (because
# the Sun's equatorial plane is tilted relative to the ecliptic plane that contains
# Earth's orbit). Disk center as seen by an observer has the same heliographic
# latitude as the heliographic latitude of the observer, and this value is known as
# the B0 angle, and varies over the year.
#
# To have no no offset (i.e., the heliographic equator is a horizontal line), you
# can specify an observer on the Sun's equatorial plane instead of using Earth as
# the observer.

# sphinx_gallery_defer_figures

observer_sun = SkyCoord(0*u.deg, 0*u.deg, 1*u.AU, obstime='2013-10-28',
                    frame=frames.HeliographicStonyhurst)
skycoord_sun = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime='2013-10-28',
                    observer=observer_sun, frame=frames.Helioprojective)

################################################################################
# Place and mark coordinates on the plot.

coords = SkyCoord(xc, yc, frame=blank_map.coordinate_frame)
ax.plot_coord(coords, 'o')
ax.plot_coord(skycoord_sun, 'o')
ax.set_title('Plotting fixed points on a blank map')

plt.show()
