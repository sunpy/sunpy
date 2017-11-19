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

#########################################################
# We first create the Map using the sample data.

aiamap = sunpy.map.Map(AIA_193_IMAGE)

#########################################################
# Now we do a quick plot

plt.figure()
aiamap.plot()
plt.colorbar()
plt.show()


####################################################################
# We now plot the aiamap data to see the distribution of the data.

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111, projection='3d')
x = y = np.arange(0, len(aiamap.data))
X, Y = np.meshgrid(x, y)
zs = np.array([aiamap.data[x][y] for x,y in zip(np.ravel(X), np.ravel(Y))])
Intensity = zs.reshape(X.shape)
ax.plot_surface(X, Y, Intensity)
ax.set_xlabel('X Coordinates')
ax.set_ylabel('Y Coordinates')
ax.set_zlabel('Intensity')
plt.show()




#####################################################################################
# Next we check for the position in the aiamap data where we get local maxima.
# We first set a threshold value, such that we will consider only those peaks greater 
# than the threshold value.

thres = 0.3
thres = thres * (np.max(Intensity) - np.min(Intensity)) + np.min(Intensity)
indexes = argrelmax(Intensity)
filter_ind = ind = np.where(Intensity[indexes] > thres )


############################################################################
# We now store the x , y coordinates where we get such local peaks

xmax = indexes[0][ind[0]]
ymax= indexes[1][ind[0]]


############################################################################
# We now check for the indices at which we get such a local maxima and plot
# those positions marked red in the aiamap data.

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111, projection='3d')
zs2 = np.array([Intensity[x][y] for x,y in zip(xmax, ymax)])
ax.plot_surface(X, Y, Intensity)
ax.scatter(ymax,xmax,zs2,color = 'r')
ax.set_xlabel('X Coordinates')
ax.set_ylabel('Y Coordinates')
ax.set_zlabel('Intensity')
plt.show()


###################################################################################
# We therefore import the coordinate functionality.

hpc_max = aiamap.pixel_to_data(xmax*u.pixel, ymax*u.pixel)


################################################################################
# Finally we do an aia plot to check for the specified locations (marked x)
# corresponding to local peaks.

fig = plt.figure()
ax = plt.subplot(projection=aiamap)
aiamap.plot()
ax.plot_coord(hpc_max, 'bx')
plt.show()

