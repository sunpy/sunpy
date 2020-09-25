#!/usr/bin/env python
# coding: utf-8

# # Using a Custom Affine Transform in the Sunpy Map rotate() function
# #### (with OpenCV)
# ----------
# How to construct a custom affine transform for in `map.rotate`.
# Requires OpenCV (cv2) Python library.

import sunpy.data.sample
import sunpy.map
import numbers
import numpy as np
import cv2

#############################
# Rotating a map in sunpy (via `map.rotate()`) uses the `skimage` or `scipy` libraries by default.
# However, the `method=` argument in `sunpy.map.rotate` can accept a custom function designed
# to use an external library. Here, we will define an affine transform function
# with the OpenCV library.
# Our custom method must have a similar function call to `sunpy.image.transform.affine_transform`
# and return the same output, namely the rotated, scaled, and transformed image.
# The required arguments are:

# Parameters
# ----------
# image: `numpy.ndarray`
#     2D image to be rotated
# rmatrix : `numpy.ndarray` that is 2x2
#     Linear transformation rotation matrix.
# order : `int` 0-5
#     Interpolation order to be used
# scale : `float`
#     A scale factor for the image
# missing : `float`
#    The value to replace any missing data after the transformation.
# image_center : tuple, optional
#     The point in the image to rotate around (axis of rotation).
#     Defaults to the center of the array.
# recenter : `bool` or array-like, optional
#     Move the axis of rotation to the center of the array or recenter coords.
#     Defaults to `True` i.e., recenter to the center of the array.


# NOTE: the required libraries have already been imported globally
# (see above: cv2, np, numbers)
def cv_rotate(image, rmatrix, order, scale, missing, image_center, recenter):
    """
    Uses `cv2.warpAffine` to do the affine transform on input `image` in same manner
    as sunpy's default `skimage.transform.warp`.
    """

    # Flags for converting input order from `integer` to the appropriate interpolation flag
    # As of Sept. 2020, OpenCV warpAffine does not support order 2,4,5
    _CV_ORDER_FLAGS = {
        0: cv2.INTER_NEAREST,
        1: cv2.INTER_LINEAR,
        2: cv2.INTER_CUBIC,
        3: cv2.INTER_CUBIC,
        4: cv2.INTER_CUBIC,
        5: cv2.INTER_CUBIC
    }

    # convert order to appropriate cv2 flag
    order = _CV_ORDER_FLAGS[max(0, min(order, 5))]

    # needed to convert `missing` from potentially a np.dtype
    # to the native `int` type required for cv2.warpAffine
    try:
        missing = missing.tolist()
    except AttributeError:
        pass

    # calculate image shift in the same manner as skimage/scipy methods
    # adapted from sunpy.transform source code

    # Make sure the image center is an array and is where it's supposed to be
    array_center = (np.array(image.shape)[::-1] - 1) / 2.0

    if image_center is not None:
        image_center = np.asanyarray(image_center)
    else:
        image_center = array_center

    # Determine center of rotation based on use (or not) of the recenter keyword
    if recenter:
        rot_center = array_center
    else:
        rot_center = image_center

    # OpenCV applies the shift+rotation operations in a different order(?); we need to calculate
    # translation using `rmatrix/scale`, but scale+rotation with `rmatrix*scale`
    # in order to match what skimage/scipy do
    displacement = np.dot(rmatrix/scale, rot_center)
    shift = image_center - displacement
    ###

    # get appropriate cv transform matrix
    # (with a slight amount of voodoo to adjust for different coordinate systems)
    rmatrix = rmatrix*scale

    trans = np.eye(3, 3)
    rot_scale = np.eye(3, 3)
    trans[:2, 2] = [-shift[0], -shift[1]]
    rot_scale[:2, :2] = rmatrix.T
    rmatrix = (rot_scale @ trans)[:2]

    # cast input image to float, if needed
    # code adapted from sunpy.transform source code
    if issubclass(image.dtype.type, numbers.Integral):
        adjusted_image = image.astype(np.float64)
    else:
        adjusted_image = image.copy()

    h, w = adjusted_image.shape

    # equivalent to skimage.transform.warp(adjusted_image, tform, order=order,
    #                                     mode='constant', cval=adjusted_missing)
    return cv2.warpAffine(adjusted_image, rmatrix, (w, h), flags=order,
                          borderMode=cv2.BORDER_CONSTANT, borderValue=missing)


# Now we want to test our implementation!
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

map_r_sk = aia_map.rotate(order=3, recenter=True, method='skimage')
map_r_sci = aia_map.rotate(order=3, recenter=True, method='scipy')
map_r_cv = aia_map.rotate(order=3, recenter=True, method=cv_rotate)


# Since skimage, scipy, and cv all use different algorithms, we need a basis for comparison:
# what is an acceptable margin of difference in the data arrays?
# We can use the Symmetric Mean Absolute Percentage Error
# (https://en.wikipedia.org/wiki/Symmetric_mean_absolute_percentage_error;
# the third equation in the link). This is a bit rudimentary, since we are looking
# at all of the pixels (and not just the solar disk), but it works for a quick approach.

# SymmetricMeanAbsolutePercentageError
# returns between 0-100% of average error between arr1 and arr2
def smape(arr1, arr2):
    eps = np.finfo(np.float64).eps
    return ((np.abs(arr1-arr2) / (np.maximum(np.abs(arr1)+np.abs(arr2), eps))).mean() * 100)


print("sk vs. sci", smape(map_r_sk.data, map_r_sci.data))
print("sci vs. cv", smape(map_r_sci.data, map_r_cv.data))
print("sk vs. cv", smape(map_r_sk.data, map_r_cv.data))

# This essentially means that the average difference between every pixel in the `skimage` and
# `scipy` images is about 1.9\%, and the `cv2` implementation results in about a 1.3-1.4\%
# average difference in the final rotated image. So, our `openCV` implementation is at least within
# the accepted difference between `scipy` and `skimage`.

# A visual comparison:

map_r_sk.peek()
map_r_sci.peek()
map_r_cv.peek()
