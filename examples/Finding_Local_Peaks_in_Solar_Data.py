"""
=================================
Finding Local Peaks in Solar Data
=================================

How to find regions of local maxima
"""
import astropy.units as u
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import sunpy.map
from sunpy.data.sample import AIA_193_IMAGE

###############################################################################
# We will first create a Map using some sample data and display it.

aiamap = sunpy.map.Map(AIA_193_IMAGE)
plt.figure()
aiamap.plot()
plt.colorbar()


###############################################################################
# Before we find regions of local maxima, we need to create some variables that
# store pixel coordinates for the 2D SDO/AIA data we are using.
# These variables are used for plotting in 3D later on.

x = y = np.arange(0, len(aiamap.data))
X, Y = np.meshgrid(x, y)

###############################################################################
# We will only consider peaks within the AIA data that are above a threshold
# value. Once these values are thresholded, the next step is to
# calcualte the pixel locations of local maxima posistions.

thres = 0.3 * (np.max(aiamap.data) - np.min(aiamap.data)) + np.min(aiamap.data)
idxs = argrelmax(aiamap.data)
thres_idxs = np.where(aiamap.data[idxs] > thres)


###############################################################################
# We now store the (x, y) coordinates of the regions where we get
# the local maxima.

xmax = idxs[0][thres_idxs[0]]
ymax = idxs[1][thres_idxs[0]]


###############################################################################
# We now check for the indices at which we get such a local maxima and plot
# those positions marked red in the aiamap data.

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, aiamap.data)
ax.view_init(elev=39, azim=64)
zs2 = np.array([aiamap.data[x][y] for x, y in zip(xmax, ymax)])
ax.scatter(ymax, xmax, zs2, color='r')
ax.set_xlabel('x coordinate')
ax.set_ylabel('y coordinate')
ax.set_zlabel('Intensity')


###############################################################################
# Now we need to turn the pixel coordinates into the world location so
# they can be easily overlaid on the Map.

hpc_max = aiamap.pixel_to_world(ymax*u.pixel, xmax*u.pixel)


###############################################################################
# Finally we do an AIA plot to check for the local maxima locations
# which will be marked with a `x` label.

fig = plt.figure()
ax = plt.subplot(projection=aiamap)
aiamap.plot()
ax.plot_coord(hpc_max, 'bx')
plt.show()
