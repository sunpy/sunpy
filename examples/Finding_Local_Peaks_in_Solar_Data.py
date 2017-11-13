"""
=================================
Finding Local Peaks in Solar Data
=================================

How to find regions of local maxima
"""
import astropy.units as u

from scipy.signal import argrelmax
import matplotlib.pyplot as plt
import numpy as np

import sunpy.map
from sunpy.data.sample import AIA_193_IMAGE

#########################################################
# We first create the Map using the sample data and then 

aiamap = sunpy.map.Map(AIA_193_IMAGE)

#########################################################
# Now we do a plot

plt.figure()
aiamap.plot()
plt.colorbar()
plt.show()


####################################################################
# We now plot the aiamap data to see the distribution of the data.
# We store the data values in yvalues.

yvalues = (aiamap.data).reshape(-1)
xvalues = [i for i in range(len(yvalues))]

# We simply plot the data

plt.figure(figsize=(10,6))
plt.plot(xvalues, yvalues)
plt.show()


#####################################################################################
# Next we check for the position in the aiamap data where we get local maxima.
# We first set a threshold value, such that we will consider only those peaks greater 
# than the threshold value.

thres = 0.3
thres = thres * (np.max(yvalues) - np.min(yvalues)) + np.min(yvalues)
indexes = argrelmax(yvalues)
filter_ind = np.where(yvalues[indexes] > thres )


############################################################################
# We now check for the indices at which we get such a local maxima and plot
# those positions marked 'x' in the aiamap data.

plt.figure(figsize=(10,6))
plt.plot(np.array(xvalues), yvalues)
plt.plot(indexes[0][filter_ind[0]] ,yvalues[indexes[0][filter_ind[0]]],'bx',c='r')
plt.title('First estimate')
plt.show()


###################################################################################
# We therefore import the coordinate functionality.

max_indices = np.unravel_index(indexes[0][filter_ind[0]], aiamap.data.shape) * u.pixel
hpc_max = aiamap.pixel_to_data(max_indices[1], max_indices[0])


################################################################################
# Finally we do an aia plot to check for the specified locations (marked x)
# corresponding to local peaks.

fig = plt.figure()
ax = plt.subplot(projection=aiamap)
aiamap.plot()
ax.plot_coord(hpc_max, 'bx')
plt.show()
