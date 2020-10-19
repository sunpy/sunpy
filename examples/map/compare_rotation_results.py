# # Comparing Results of `map.rotate` Between Library Implementations
# ----------

import sunpy.data.sample
import sunpy.map
import numbers
import numpy as np
import cv2

# Rotating a map in Sunpy has a choice between three libraries: `skimage', `scipy`, and `cv2`.
# One can also create a custom rotation function and drop it into `map.rotate`.
# Since all these options use different algorithms, we need a basis for comparison:
# what is an acceptable margin of difference in the final, rotated data product?
# One option is to use the Symmetric Mean Absolute Percentage Error
# (https://en.wikipedia.org/wiki/Symmetric_mean_absolute_percentage_error;
# the third equation in the link). This is a bit rudimentary, since it looks
# at all of the pixels in the image (and not just the solar disk), but it works
# for a quick check.

# SymmetricMeanAbsolutePercentageError
# returns between 0%-100% of average error between arr1 and arr2
def smape(arr1, arr2):
    eps = np.finfo(np.float64).eps
    return ((np.abs(arr1-arr2) / (np.maximum(np.abs(arr1)+np.abs(arr2), eps))).mean() * 100)

# get sample data, and do a default rotation based on FITS header data.
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

map_r_sk = aia_map.rotate(order=3, recenter=True, method='skimage')
map_r_sci = aia_map.rotate(order=3, recenter=True, method='scipy')
map_r_cv = aia_map.rotate(order=3, recenter=True, method='cv2')

# Calculate Symmetric Mean Percentage Error
print("sk vs. sci", smape(map_r_sk.data, map_r_sci.data))
print("sci vs. cv", smape(map_r_sci.data, map_r_cv.data))
print("sk vs. cv", smape(map_r_sk.data, map_r_cv.data))

# This essentially means that the average difference between every pixel in the `skimage` and
# `scipy` images is about 1.9\%, and the `cv2` implementation results in about a 1.3-1.4\%
# average difference in the final rotated image. So, the `openCV` implementation is at least within
# the accepted difference between `scipy` and `skimage`.

# A visual comparison:

map_r_sk.peek()
map_r_sci.peek()
map_r_cv.peek()
