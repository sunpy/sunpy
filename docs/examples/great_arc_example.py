"""
============================
Using a Great Arc
============================
How to draw a great arc and extract the data along the arc
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
# Let's first consider a simple example of transforming a single coordinate.
# Here we consider the "center" coordinate of 0,0 degrees
start = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=m.coordinate_frame)
end = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=m.coordinate_frame)

###############################################################################
# Create the great arc between the start and end points
great_arc = GreatArc(start, end)

###############################################################################
# Plot the great arc on the Sun
fig = plt.figure()
ax = plt.subplot(projection=m)
m.plot(axes=ax)
ax.plot_coord(great_arc.coordinates(), color='c')
plt.show()

###############################################################################
# Now get the pixels out
pixels = np.asarray(np.rint(great_arc.coordinates().to_pixel(m.wcs)), dtype=int)
fig, ax = plt.subplots()
ax.plot(great_arc.inner_angles().to(u.deg), m.data[pixels[0, :], pixels[1, :]])
ax.set_xlabel('degrees of arc from initial position')
ax.set_ylabel('intensity')
