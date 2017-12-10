"""
=================================
Finding Local Peaks in Solar Data
=================================

How to find regions of local maxima
"""
import astropy.units as u
from skimage.feature import peak_local_max
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
# This function comes from sci-kit image and the documenation is found
# here `~skimage.feature.peak_local_max`.

coordinates = peak_local_max(aiamap.data, min_distance=60, threshold_rel=0.2)


###############################################################################
# We now check for the indices at which we get such a local maxima and plot
# those positions marked red in the aiamap data.

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, aiamap.data)
ax.view_init(elev=39, azim=64)
zs2 = np.array([aiamap.data[x][y] for x, y in zip(coordinates[:, 0], coordinates[:, 1])])
ax.scatter(coordinates[:, 1], coordinates[:, 0],zs2, color = 'r')
ax.set_xlabel('X Coordinates')
ax.set_ylabel('Y Coordinates')
ax.set_zlabel('Intensity')


###############################################################################
# Now we need to turn the pixel coordinates into the world location so
# they can be easily overlaid on the Map.

hpc_max = aiamap.pixel_to_data(coordinates[:, 1]*u.pixel, coordinates[:, 0]*u.pixel)


###############################################################################
# Finally we do an AIA plot to check for the local maxima locations
# which will be marked with a `x` label.

fig = plt.figure()
ax = plt.subplot(projection=aiamap)
aiamap.plot()
ax.plot_coord(hpc_max, 'bx')
plt.show()
