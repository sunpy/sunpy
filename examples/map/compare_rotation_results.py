"""
================================
Comparing Map Rotation Functions
================================

This example will compare between the current library implementations for `sunpy.map.GenericMap.rotate`.
"""
import matplotlib.colors
import numpy as np
import pylab as plt

import astropy.units as u

import sunpy.data.sample
import sunpy.map

###############################################################################
# Rotating a map in sunpy has a choice between three libraries: ``scipy`` (the default),
# ``skimage`` and ``opencv2``. Furthermore, one can also create a custom rotation
#  function and register it for use with `~sunpy.map.GenericMap.rotate`, see :ref:`map_rotate_custom`.
#
# Since all these options use different algorithms, we need a basis for comparison.
#
# The problem is, defining a acceptable margin of difference in the final, rotated data product
# is difficult. One option is to use the `Symmetric Mean Absolute Percentage Difference
# <https://en.wikipedia.org/wiki/Symmetric_mean_absolute_percentage_error>`__
#
# This is a bit rudimentary, since it looks at all of the pixels in the image
# (and not just the solar disk), but it works as a quick check.


def smape(arr1, arr2):
    """
    Symmetric Mean Absolute Percentage Error.

    Parameters
    ----------
    arr1 : `numpy.ndarray`
        First array to compare.
    arr2 : `numpy.ndarray`
        Second array to compare.

    Returns
    -------
    `float`
        Between 0% - 100% of average error between ``arr1`` and ``arr2``
    """
    eps = np.finfo(np.float64).eps
    return np.nanmean(np.abs(arr1 - arr2) / (np.maximum(np.abs(arr1) + np.abs(arr2), eps))) * 100


###############################################################################
# Using an HMI sample data, we will do a rotation to align the image to the north.
# By default, the order of rotation is 3.

hmi_map = sunpy.map.Map(sunpy.data.sample.HMI_LOS_IMAGE)

scipy_map = hmi_map.rotate(method='scipy')
skimage_map = hmi_map.rotate(method='scikit-image')
cv2_map = hmi_map.rotate(method='opencv')

###############################################################################
# Now we will calculate the Symmetric Mean Percentage Error.

print("HMI Example")
print(f"scipy vs scikit-image {smape(scipy_map.data, skimage_map.data):.2f}%")
print(f"scipy vs opencv2 {smape(scipy_map.data, cv2_map.data):.2f}%")
print(f"scikit-image vs opencv2 {smape(skimage_map.data, cv2_map.data):.2f}%")

###############################################################################
# This essentially means that the average difference between every pixel in the
# ``scipy`` and ``skimage`` rotated images are about 15\%.
# Where as the comparision with ``opencv2`` is only about 10\%.

###############################################################################
# Now for a visual comparison, the raw differences, that should highlight the differences.
# Note that only two comparisions are shown.

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
norm = matplotlib.colors.SymLogNorm(5, vmin=-30, vmax=30)

img1 = ax1.imshow(scipy_map.data - skimage_map.data, cmap='RdBu_r', norm=norm)
ax1.set_title("HMI Difference: scipy vs scikit-image")
fig.colorbar(img1, ax=ax1)

img2 = ax2.imshow(scipy_map.data - cv2_map.data, cmap='RdBu_r', norm=norm)
ax2.set_title("HMI Difference: scipy vs opencv2")
fig.colorbar(img2, ax=ax2)

plt.show()

###############################################################################
# We can repeat this but for AIA data, using a 171 sample image.
# Since there is no rotation in the FITS header, we will rotate by 30 degrees

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

scipy_map = aia_map.rotate(30*u.deg, method='scipy')
skimage_map = aia_map.rotate(30*u.deg, method='scikit-image')
cv2_map = aia_map.rotate(30*u.deg, method='opencv')

###############################################################################
# Now we will calculate the Symmetric Mean Percentage Error for AIA.

print("AIA 171 Example")
print(f"scipy vs scikit-image {smape(scipy_map.data, skimage_map.data):.2f}%")
print(f"scipy vs opencv2 {smape(scipy_map.data, cv2_map.data):.2f}%")
print(f"scikit-image vs opencv2 {smape(skimage_map.data, cv2_map.data):.2f}%")

###############################################################################
# This essentially means that the average difference between every pixel in the
# ``scipy`` and ``skimage`` rotated images are about 2\%.
# Where as the comparision with ``opencv2`` is only about 1.5\%.

###############################################################################
# Now for a visual comparison, the raw differences, that should highlight the differences.
# Note that only two comparisions are shown.

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
norm = matplotlib.colors.SymLogNorm(5, vmin=-30, vmax=30)

img1 = ax1.imshow(scipy_map.data - skimage_map.data, cmap='RdBu_r', norm=norm)
ax1.set_title("AIA Difference: scipy vs scikit-image")
fig.colorbar(img1, ax=ax1)

img2 = ax2.imshow(scipy_map.data - cv2_map.data, cmap='RdBu_r', norm=norm)
ax2.set_title("AIA Difference: scipy vs opencv2")
fig.colorbar(img2, ax=ax2)

plt.show()
