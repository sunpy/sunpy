"""
============================
Longitude and Latitude Lines
============================

How to draw your own (Stonyhurst) longitude and latitude lines
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from sunpy.data.sample import AIA_171_IMAGE
from sunpy.map import Map

###############################################################################
# We first create the Map using the sample data and import the coordinate
# functionality.
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
aia = Map(AIA_171_IMAGE)

###############################################################################
# Let's first consider a simple example of transforming a single coordinate.
# Here we consider the "center" coordinate of 0,0 degrees
stonyhurst_center = SkyCoord(0 * u.deg, 0 * u.deg, frame=frames.HeliographicStonyhurst)

###############################################################################
# Next we transform it into the coordinate frame of our map which is in
# helioprojective coordinates though we don't really need to know that.
hpc_stonyhurst_center = stonyhurst_center.transform_to(aia.coordinate_frame)
print(hpc_stonyhurst_center)

###############################################################################
# Now let's consider transform two lines, one of zero longitude and one of
# of zero latitude. First define the coordinates as we did before and then
# transform them.
lon0 = SkyCoord(np.linspace(-80, 80, 100) * u.deg,
                np.zeros(100) * u.deg, frame=frames.HeliographicStonyhurst)
lat0 = SkyCoord(np.zeros(100)*u.deg,
                np.linspace(-90, 90, 100)*u.deg, frame=frames.HeliographicStonyhurst)

hpc_lon0 = lon0.transform_to(aia.coordinate_frame)
hpc_lat0 = lat0.transform_to(aia.coordinate_frame)

###############################################################################
# Let's now plot the results
fig = plt.figure()
ax = fig.add_subplot(111)
aia.plot(axes=ax)
ax.plot(hpc_lon0.Tx, hpc_lon0.Ty, color='white')
ax.plot(hpc_lat0.Tx, hpc_lat0.Ty, color='white')
plt.savefig("lonlat0.png")
#plt.show()
