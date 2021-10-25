"""
================================
Comparing Map Rotation Functions
================================

Comparing between library implementations for `map.rotate`.
"""

import matplotlib.pyplot as plt
import numpy as np

import sunpy.data.sample
import sunpy.map


###############################################################################
# Rotating a map in sunpy has a choice between two libraries: `skimage' and `scipy`
# One can also create a custom rotation function and drop it into `map.rotate`.
# Since all these options use different algorithms, we need a basis for comparison:
# what is an acceptable margin of difference in the final, rotated data product?
# One option is to use the Symmetric Mean Absolute Percentage Difference
# (https://en.wikipedia.org/wiki/Symmetric_mean_absolute_percentage_error;
# the third equation in the link). This is a bit rudimentary, since it looks
# at all of the pixels in the image (and not just the solar disk), but it works
# for a quick check.
def smape(arr1, arr2):
    """
    SymmetricMeanAbsolutePercentageError
    returns between 0%-100% of average error between arr1 and arr2
    """
    eps = np.finfo(np.float64).eps
    return ((np.abs(arr1 - arr2) / (np.maximum(np.abs(arr1) + np.abs(arr2), eps))).mean() * 100)


###############################################################################
# Get sample data, and do a default rotation based on FITS header data to align
# solar North with the y-axis.
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

map_r_sk = aia_map.rotate(order=3, recenter=True, method='skimage')
map_r_sci = aia_map.rotate(order=3, recenter=True, method='scipy')

###############################################################################
# Calculate Symmetric Mean Percentage Error
print("sk vs. sci", smape(map_r_sk.data, map_r_sci.data))

# This essentially means that the average difference between every pixel in the ``skimage`` and
# ``scipy`` images is about 1.9\%.

###############################################################################
# Make some plots to compare the rotated images.
# A visual comparison:

map_r_sk.peek()
map_r_sci.peek()

# Comparing raw differences:

fig, ax = plt.subplots()
img = ax.imshow(map_r_sk.data - map_r_sci.data, vmin=-150, vmax=150, cmap=plt.cm.seismic)
fig.colorbar(img)
ax.set_title("Raw Difference: Sklearn vs. Scipy")
