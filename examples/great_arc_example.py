"""
=============================
Drawing and using a Great Arc
=============================

This example shows you how to define and draw a great arc on an image of the
Sun, and to extract intensity values along that arc from the data.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunpy.coordinates.utils import GreatArc
import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE
m = sunpy.map.Map(AIA_171_IMAGE)

###############################################################################
# Let's define the start and end co-ordinates of the arc on the Sun.
start = SkyCoord(735*u.arcsec, -471*u.arcsec, frame=m.coordinate_frame)
end = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=m.coordinate_frame)

###############################################################################
# Create the great arc between the start and end points.
great_arc = GreatArc(start, end)

###############################################################################
# Plot the great arc on the Sun
fig = plt.figure()
ax = plt.subplot(projection=m)
m.plot(axes=ax)
ax.plot_coord(great_arc.coordinates(), color='c')
plt.show()

###############################################################################
# Now we can calculate the nearest integer pixels of the data that correspond
# to the location of arc.
pixels = np.asarray(np.rint(great_arc.coordinates().to_pixel(m.wcs)), dtype=int)
x = pixels[0, :]
y = pixels[1, :]

###############################################################################
# Get the intensity along the arc from the start to the end point
intensity_along_arc = m.data[y, x]

###############################################################################
# Define the angular location of each pixel along the arc from the start point
# to the end
angles = great_arc.inner_angles().to(u.deg)

###############################################################################
#
fig, ax = plt.subplots()
ax.plot(angles, intensity_along_arc)
ax.set_xlabel('degrees of arc from start')
ax.set_ylabel('intensity')
ax.grid(linestyle='dotted')
plt.show()
