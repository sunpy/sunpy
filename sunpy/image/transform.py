"""
Functions for geometrical image transformation and warping.
"""

from __future__ import absolute_import

import warnings

import numpy as np
import scipy.ndimage.interpolation
try:
    import skimage.transform
    scikit_image_not_found = False
except ImportError:  # pragma: no cover
    scikit_image_not_found = True  # pragma: no cover

__all__ = ['affine_transform']

def affine_transform(image, rmatrix, order=3, scale=1.0, image_center=None,
                     recenter=False, missing=0.0, use_scipy=False):
    """
    Rotates, shifts and scales an image using :func:`skimage.transform.warp`, or
    :func:`scipy.ndimage.interpolation.affine_transform` if specified. Falls back to
    the scipy function if scikit-image can't be imported.

    Parameters
    ----------
    image : ndarray
        2D Image to be rotated.
    rmatrix : 2x2
        Linear transformation rotation matrix.
    order : int 0-5
        Interpolation order to be used. When using scikit-image this parameter
        is passed into :func:`skimage.transform.warp` (e.g., 3 corresponds to
        bi-cubic interpolation).
        When using scipy it is passed into
        :func:`scipy.ndimage.interpolation.affine_transform` where it controls
        the order of the spline.
        Default: 3
    scale : float
        A scale factor for the image. Default is no scaling.
    image_center : tuple
        The point in the image to rotate around (axis of rotation).
        Default: center of the array.
    recenter : bool or array-like
        Move the axis of rotation to the center of the array or recenter coords.
        Default: True, recenter to the center of the array.
    missing : float
        The value to replace any missing data after the transformation.
    use_scipy : bool
        Force use of :func:`scipy.ndimage.interpolation.affine_transform`.
        Will set all NaNs in image to zero before doing the transform.
        Default: False unless sckit-image is not installed.

    Returns
    -------
    out : New rotated, scaled and translated image.

    Notes
    -----
    This function will cast any integer input data to float64, it is possible to
    manually undo this by casting back with `numpy.ndarray.astype`, however to avoid
    possible loss of precision this is not done automatically.

    This function is used throughout the SunPy code base as an equivalent to
    the IDL's rot() function. However, this function does not use the same
    algorithm as the IDL rot() function.
    IDL's rot() calls the `POLY_2D <http://www.exelisvis.com/docs/poly_2d.html>`_
    method to calculate the inverse mapping of original to target pixel
    coordinates. This is a polynominal geometrical transformation.
    Then optionally it uses a bicubic convolution interpolation
    algorithm to map the original to target pixel values.

    This algorithm uses an affine transformation as opposed to a polynomial
    geometrical transformation. It then defaults to a bicubic convolution
    interpolation.

    The two affine transform routines this method,
    :func:`scipy.ndimage.interpolation.affine_transform` and
    :func:`skimage.transform.warp` do not use the same interpolation
    algorithm and therefore should not be expected to give identical output,
    the scikit-image routines are preferred.

    No direct comparision between this function and the IDL rot() function has
    been performed as it is a different transformation, however both are
    assumed to be mathamatically vaild.
    """

    rmatrix = rmatrix / scale
    array_center = (np.array(image.shape)[::-1]-1)/2.0

    # Make sure the image center is an array and is where it's supposed to be
    if image_center is not None:
        image_center = np.asanyarray(image_center)
    else:
        image_center = array_center

    # Determine center of rotation based on use (or not) of the recenter keyword
    if recenter:
        rot_center = array_center
    else:
        rot_center = image_center

    displacement = np.dot(rmatrix, rot_center)
    shift = image_center - displacement

    if use_scipy or scikit_image_not_found:
        if np.any(np.isnan(image)):
            warnings.warn("Setting NaNs to 0 for SciPy rotation", Warning)
        # Transform the image using the scipy affine transform
        rotated_image = scipy.ndimage.interpolation.affine_transform(
                np.nan_to_num(image).T, rmatrix, offset=shift, order=order,
                mode='constant', cval=missing).T
    else:
        # Make the rotation matrix 3x3 to include translation of the image
        skmatrix = np.zeros((3, 3))
        skmatrix[:2, :2] = rmatrix
        skmatrix[2, 2] = 1.0
        skmatrix[:2, 2] = shift
        tform = skimage.transform.AffineTransform(skmatrix)

        # Transform the image using the skimage function
        # Image data is normalised because warp() requires an array of values
        # between -1 and 1.
        if np.issubdtype(image.dtype, np.integer):
            adjusted_image = image.astype(np.float64)
        else:
            adjusted_image = image.copy()
        if np.any(np.isnan(adjusted_image)) and order >= 4:
            warnings.warn("Setting NaNs to 0 for higher-order scikit-image rotation", Warning)
            adjusted_image = np.nan_to_num(adjusted_image)

        im_min = np.nanmin(adjusted_image)
        adjusted_image -= im_min
        im_max = np.nanmax(adjusted_image)
        adjusted_image /= im_max
        adjusted_missing = (missing - im_min) / im_max
        rotated_image = skimage.transform.warp(adjusted_image, tform, order=order,
                                               mode='constant', cval=adjusted_missing)

        rotated_image *= im_max
        rotated_image += im_min

        if rotated_image.dtype != image.dtype:
            warnings.warn("Rotating this image has change it's dtype", RuntimeWarning)
    return rotated_image
