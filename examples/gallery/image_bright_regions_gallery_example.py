
# coding: utf-8

"""
====================================================
Manipulating Map image data - finding bright regions
====================================================

This example shows how you can do basic image processing on SunPy map image data. In this example, we try to find the brightest regions in an AIA image, and count the approximate number of regions of interest.
"""

##############################################################################
# First, import the modules we will need:


import sunpy
from sunpy import map
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

##############################################################################
# Now, we create a SunPy Map object from an AIA FITS file.


aiamap = map.Map('AIA20110607_062209_0193.fits')

##############################################################################
# Let's plot the map. Here we use the `clim` method to clip the colour table. This show features in an image with high dynamic range more clearly.


aiamap.plot()
plt.clim([0,2500])
plt.colorbar()
plt.show()

##############################################################################
# Now we want to find the brightest regions in this image. We start by finding the maximum value in the image data.


m = np.max(aiamap.data)

##############################################################################
# Now we want to make a mask, which tells us which regions are bright. We choose the criterion that the data should be at least 5% of the maximum value. Pixels with intensity # values greater than this are included in the mask, while all other pixels are excluded. 


mask = aiamap.data > m*0.05
mask

##############################################################################
# Mask is a `boolean` array. It can be used to modify the original map data.



data2 = aiamap.data * mask

##############################################################################
# We now have a new set of data, where all values outside the mask are set to zero. We can create a new SunPy map structure with this data.


aiamap2 = map.Map(data2, aiamap.meta)



aiamap2.plot()
plt.colorbar()
plt.show()

##############################################################################
# Only the brightest pixels remain in the image. However, these areas are artificially broken up into small regions. Estimating the number of significant hot regions will be
# difficult. We can solve this by applying some smoothing to the image data. Here we apply a 2D Gaussian smoothing function to the data.


data3 = ndimage.gaussian_filter(data2,32)

##############################################################################
# Now we can make a third SunPy map with this smoothed data.


aiamap3 = sunpy.map.Map(data3, aiamap.meta)

##############################################################################
# Now we can use the function `label` from the `scipy.ndimage` module, which counts the number of contiguous regions in an image.


labels, n = ndimage.label(aiamap3.data)

##############################################################################
# Finally, we plot the smoothed bright image data, along with the estimate of the number of distinct regions. We can see that approximately 6 distinct hot regions are present # above the 5% of the maximum level.


aiamap3.plot()
plt.figtext(0.3,0.2,'Number of regions = ' + str(n),color='white')
plt.show()

