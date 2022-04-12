"""
================================
Comparing Map Rotation Functions
================================

This example will compare between the current library implementations for `sunpy.map.GenericMap.rotate`.
"""
import numpy as np
import pylab as plt

import sunpy.data.sample
import sunpy.map

###############################################################################
# Rotating a map in sunpy has a choice between three libraries: `skimage', `scipy`,
# and `cv2`. One can also create a custom rotation function and add it, see :ref:`map_rotate_custom`.
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
    Symmetric Mean Absolute Percentage Error

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
# We will also fix the order (to 3) for each operation.


hmi_map = sunpy.map.Map(sunpy.data.sample.HMI_LOS_IMAGE)

skimage_map = hmi_map.rotate(order=3, recenter=True, method='skimage')
scipy_map = hmi_map.rotate(order=3, recenter=True, method='scipy')
cv2_map = hmi_map.rotate(order=3, recenter=True, method='cv2')

###############################################################################
# Now we will calculate the Symmetric Mean Percentage Error.

print("HMI Example")
print(f"scikit-image vs scipy {smape(skimage_map.data, scipy_map.data):.2f}%")
print(f"scipy vs OpenCV {smape(scipy_map.data, cv2_map.data):.2f}%")
print(f"scikit-image vs OpenCV {smape(skimage_map.data, cv2_map.data):.2f}%")

###############################################################################
# This essentially means that the average difference between every pixel in the ``skimage`` and
# ``scipy`` images is about 15\%, and the ``cv2`` implementation results in about a 100\%
# average difference in the final rotated image. In this case ``openCV`` implementation
# can not handle the HMI data very well.

###############################################################################
# Now for a visual comparison, that should highlight the differences.

skimage_map.plot()
plt.show()

scipy_map.plot()
plt.show()

cv2_map.plot()
plt.show()

###############################################################################
# Furthermore, we can do the raw differences between the images.

fig0, ax0 = plt.subplots()
img0 = ax0.imshow(skimage_map.data - scipy_map.data, vmin=-150, vmax=150, cmap=plt.cm.seismic)
fig0.colorbar(img0)
ax0.set_title("Raw Difference: scikit-image vs scipy")

plt.show()

fig1, ax1 = plt.subplots()
img1 = ax1.imshow(scipy_map.data - cv2_map.data, vmin=-150, vmax=150, cmap=plt.cm.seismic)
fig1.colorbar(img1)
ax1.set_title("Raw Difference: scipy vs OpenCV")

plt.show()

fig2, ax2 = plt.subplots()
img2 = ax2.imshow(skimage_map.data - cv2_map.data, vmin=-150, vmax=150, cmap=plt.cm.seismic)
fig2.colorbar(img2)
ax2.set_title("Raw Difference: scikit-image vs OpenCV")

plt.show()

###############################################################################
# We can repeat this but for AIA data, using a 171 sample image.
# Since there is no rotation in the FITS header, it will just apply the same rotation to the image.

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

skimage_map = aia_map.rotate(order=3, recenter=True, method='skimage')
scipy_map = aia_map.rotate(order=3, recenter=True, method='scipy')
cv2_map = aia_map.rotate(order=3, recenter=True, method='cv2')

###############################################################################
# Now we will calculate the Symmetric Mean Percentage Error.

print("AIA 171 Example")
print(f"scikit-image vs scipy {smape(skimage_map.data, scipy_map.data):.2f}%")
print(f"scipy vs OpenCV {smape(scipy_map.data, cv2_map.data):.2f}%")
print(f"scikit-image vs OpenCV {smape(skimage_map.data, cv2_map.data):.2f}%")

###############################################################################
# This essentially means that the average difference between every pixel in the ``skimage`` and
# ``scipy`` images is about 1.9\%, and the ``cv2`` implementation results in about a 1.3-1.4\%
# average difference in the final rotated image. So, the ``openCV`` implementation is at least within
# the accepted difference between ``scipy`` and ``skimage``.

###############################################################################
# Now for a visual comparison, that should highlight the differences.

skimage_map.plot()
plt.show()

scipy_map.plot()
plt.show()

cv2_map.plot()
plt.show()

###############################################################################
# Furthermore, we can do the raw differences between the images.

fig0, ax0 = plt.subplots()
img0 = ax0.imshow(skimage_map.data - scipy_map.data, vmin=-150, vmax=150, cmap=plt.cm.seismic)
fig0.colorbar(img0)
ax0.set_title("Raw Difference: scikit-image vs scipy")

plt.show()

fig1, ax1 = plt.subplots()
img1 = ax1.imshow(scipy_map.data - cv2_map.data, vmin=-150, vmax=150, cmap=plt.cm.seismic)
fig1.colorbar(img1)
ax1.set_title("Raw Difference: scipy vs OpenCV")

plt.show()

fig2, ax2 = plt.subplots()
img2 = ax2.imshow(skimage_map.data - cv2_map.data, vmin=-150, vmax=150, cmap=plt.cm.seismic)
fig2.colorbar(img2)
ax2.set_title("Raw Difference: scikit-image vs OpenCV")

plt.show()
