  
"""
==============================================
Plot Positions on a blank Map
==============================================
This example shows how to plot random positions on a blank map
"""
# Start by importing the necessary modules.

import matplotlib.pyplot as plt
import numpy as np

import sunpy.map
from sunpy.coordinates import frames

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

################################################################################
# Creating blank map data with zeros
data = np.zeros((1000, 1000))

# Define coordinates and frame of reference and make the header using sunpy.map.make_fitswcs_header
skycoord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime = '2013-10-28', observer = 'earth', frame = frames.Helioprojective)
header = sunpy.map.make_fitswcs_header(data, skycoord)

# Use the sunpy.map.Map to make the blank map
blank_map = sunpy.map.Map(data, header)

################################################################################

# Initialize the plot and add the map to it
fig = plt.figure()
ax = plt.subplot(projection=blank_map)
im = blank_map.plot()

# Prevent the image from being re-scaled while overplotting.
ax.set_autoscale_on(False)

# Coordinates that are being plotted - (0, 0), (50, 100) and (400, 400)
xc = [0,50,400] * u.arcsec
yc = [0,100,400] * u.arcsec

# Place and mark coordinates on the plot
coords = SkyCoord(xc, yc, frame=blank_map.coordinate_frame)
p = ax.plot_coord(coords, 'o')

# Set title
ax.set_title('Custom plot with WCSAxes')

# Add the colorbar and display the plot
plt.colorbar()
plt.show()