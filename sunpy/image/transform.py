"""
Functions for geometrical image transformation and warping.
"""
import numbers

import numpy as np
import scipy.ndimage
from scipy.signal import convolve2d

from sunpy.util.exceptions import warn_deprecated, warn_user

__all__ = ['affine_transform']


def affine_transform(image, rmatrix, order=3, scale=1.0, image_center=None,
                     recenter=False, missing=0.0, use_scipy=None, *, method='scipy', clip=True):
    """
    Rotates, shifts and scales an image.

    Parameters
    ----------
    image : `numpy.ndarray`
        2D image to be rotated.
    rmatrix : `numpy.ndarray` that is 2x2
        Linear transformation rotation matrix.
    order : `int` 0-5, optional
        Interpolation order to be used, defaults to 3.  The precise meaning depends
        on the rotation method specifed by ``method``.
    scale : `float`
        A scale factor for the image with the default being no scaling.
    image_center : tuple, optional
        The point in the image to rotate around (axis of rotation).
        Defaults to the center of the array.
    recenter : `bool` or array-like, optional
        Move the axis of rotation to the center of the array or recenter coords.
        Defaults to `True` i.e., recenter to the center of the array.
    missing : `float`, optional
        The value to replace any missing data after the transformation.
    method : {``'skimage'``, ``'scipy'``}, optional
        Transform function to use. Currently
        :func:`scipy.ndimage.affine_transform` and
        :func:`skimage.transform.warp` are supported.
        Defaults to 'scipy'.
    clip : `bool`, optional
        If `True`, clips the pixel values of the output image to the range of the
        input image (including the value of ``missing``).
        Defaults to `True`.

    Returns
    -------
    `numpy.ndarray`:
        New rotated, scaled and translated image.

    Notes
    -----
    If there are NaNs in the image, pixels in the output image will be set to NaN if
    they are within a number of pixels equal to half of the ``order`` parameter.

    For rotation using ``scikit-image``, an input image with integer data is cast to
    64-bit floats, and the output image can be re-cast using `numpy.ndarray.astype`
    if desired.
    """
    rmatrix = rmatrix / scale
    array_center = (np.array(image.shape)[::-1] - 1) / 2.0

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

    # Check for any NaNs and set them to `missing` for the rotation
    has_nans = np.any(np.isnan(image))
    adjusted_image = np.nan_to_num(image, nan=missing) if has_nans else image

    method = _get_transform_method(method, use_scipy)
    if method == 'scipy':
        # Transform using scipy, with any NaNs set to zero
        rotated_image = scipy.ndimage.affine_transform(
            adjusted_image.T, rmatrix, offset=shift, order=order,
            mode='constant', cval=missing).T

        if clip:
            # Clip the image to the input range, and assume that the `missing` value has been used
            rotated_image.clip(np.min([missing, np.min(adjusted_image)]),
                               np.max([missing, np.max(adjusted_image)]),
                               out=rotated_image)
    else:  # method == 'skimage'
        import skimage.transform

        # Make the rotation matrix 3x3 to include translation of the image
        skmatrix = np.zeros((3, 3))
        skmatrix[:2, :2] = rmatrix
        skmatrix[2, 2] = 1.0
        skmatrix[:2, 2] = shift
        tform = skimage.transform.AffineTransform(skmatrix)

        if issubclass(adjusted_image.dtype.type, numbers.Integral):
            warn_user("Integer input data has been cast to float64.")
            adjusted_image = adjusted_image.astype(np.float64)
        else:
            adjusted_image = adjusted_image.copy()

        # Scale image to range [0, 1]
        im_min = np.min(adjusted_image)
        adjusted_image -= im_min
        im_max = np.max(adjusted_image)
        adjusted_missing = missing - im_min
        if im_max > 0:
            adjusted_image /= im_max
            adjusted_missing /= im_max

        rotated_image = skimage.transform.warp(adjusted_image, tform, order=order,
                                               mode='constant', cval=adjusted_missing, clip=clip)

        # Convert the image back to its original range
        if im_max > 0:
            rotated_image *= im_max
        rotated_image += im_min

    # Restore the NaNs
    if has_nans:
        # Use a convolution to find all pixels that are affected by NaNs
        # We want a kernel size that is an odd number that is at least order+1
        size = 2*int(np.ceil(order/2))+1
        expanded_nans = convolve2d(np.isnan(image).astype(int),
                                   np.ones((size, size)).astype(int),
                                   mode='same')
        rotated_nans = scipy.ndimage.affine_transform(expanded_nans.T, rmatrix,
                                                      offset=shift, order=1,
                                                      mode='nearest').T
        rotated_image[rotated_nans > 0] = np.nan

    return rotated_image


def _get_transform_method(method, use_scipy):
    # This is re-used in affine_transform and GenericMap.rotate
    supported_methods = {'scipy', 'skimage'}
    if method not in supported_methods:
        raise ValueError(f'Method {method} not in supported methods: {supported_methods}')

    if use_scipy is not None:
        warn_deprecated("The 'use_scipy' argument is deprecated. "
                        "Specify the rotation method to the 'method' "
                        "keyword argument instead.")
        if use_scipy is True and method != 'scipy':
            warn_user(f"Using scipy instead of {method} for rotation.")
            method = 'scipy'

    if method == 'skimage':
        try:
            import skimage  # NoQA
        except ImportError:
            raise ImportError("scikit-image must be installed to be usable for rotation.")

    return method
