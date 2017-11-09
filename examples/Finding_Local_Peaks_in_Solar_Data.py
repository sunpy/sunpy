"""
=================================
Finding Local Peaks in Solar Data
=================================

How to find regions of local maxima
"""
import astropy.units as u

import scipy.signal
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
# We store the data values in y 
yvalues = (aiamap.data).reshape(-1)
xvalues = [i for i in range(len(yvalues))]

# We simply plot the data
plt.figure(figsize=(10,6))
plt.plot(xvalues, yvalues)
plt.show()


###################################################################################################
# Next we check for the position in the aiamap data where we get local maxima.
# Using the below function we take the window width of 2000 at an interval of 100 to use for calculating the CWT matrix.

indexes = scipy.signal.find_peaks_cwt(yvalues, np.arange(1,2000,100), noise_perc=0.1)


##################################################################################################
# We now check for the indices at which we get such a local maxima and plot those positions marked 'x'
# in the aiamap data
print(indexes)
print(xvalues, yvalues[indexes])
plt.figure(figsize=(10,6))
plt.plot(np.array(xvalues), yvalues)
plt.plot(indexes,yvalues[indexes],'bx',c='r')
plt.title('First estimate')
plt.show()




###################################################################################################
# We therefore import the coordinate functionality.
max_indices = np.unravel_index(indexes, aiamap.data.shape) * u.pixel
hpc_max = aiamap.pixel_to_data(max_indices[1], max_indices[0])


################################################################################
# Finally we do an aia plot to check for the specified locations (marked x)
# corresponding to local peaks.
fig = plt.figure()
ax = plt.subplot(projection=aiamap)
aiamap.plot()
ax.plot_coord(hpc_max, 'bx')
plt.show()
