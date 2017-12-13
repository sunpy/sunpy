"""
=================================
Finding Local Peaks in Solar Data
=================================

Detecting radiation peaks in the solar surface is often crucial in the study of solar flares.
This example illustrates detection of those areas where there is a spike in
solar radiation intensity. We use the peak_local_max function under the scikit library
to find those regions in the map data where the intensity values form a local maxima.
Finally we plot those peaks in the original AIA plot.
"""

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from skimage.feature import peak_local_max

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

x = np.arange(aiamap.data.shape[0])
y = np.arange(aiamap.data.shape[1])
X, Y = np.meshgrid(x, y)


#######################################################################################
# We will only consider peaks within the AIA data that have minimum intensity
# value equal to threshold_rel * max(Intensity) which is 20% of the maximum intensity
# in our case. The next step is to calculate the pixel locations of local maxima
# positions where peaks are separated by atleast min_distance = 60 pixels.
# This function comes from scikit image and the documenation is found
# here `~skimage.feature.peak_local_max`.

coordinates = peak_local_max(aiamap.data, min_distance=60, threshold_rel=0.2)


###############################################################################
# We now check for the indices at which we get such a local maxima and plot
# those positions marked red in the aiamap data.

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, aiamap.data)
ax.view_init(elev=39, azim=64)
peaks_pos = aiamap.data[coordinates[:, 0], coordinates[:, 1]]
ax.scatter(coordinates[:, 1], coordinates[:, 0], peaks_pos, color='r')
ax.set_xlabel('X Coordinates')
ax.set_ylabel('Y Coordinates')
ax.set_zlabel('Intensity')


###############################################################################
# Now we need to turn the pixel coordinates into the world location so
# they can be easily overlaid on the Map.

hpc_max = aiamap.pixel_to_world(coordinates[:, 1]*u.pixel, coordinates[:, 0]*u.pixel)

###############################################################################
# Finally we do an AIA plot to check for the local maxima locations
# which will be marked with a blue `x` label.

fig = plt.figure()
ax = plt.subplot(projection=aiamap)
aiamap.plot()
ax.plot_coord(hpc_max, 'bx')
plt.show()
