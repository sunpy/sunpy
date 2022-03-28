"""
Functions for geometrical image transformation and warping.
"""
import numbers
from functools import wraps

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
        Can be `numpy.nan`.
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

    # While `use_scipy` is still supported, we have to check which method to actually use
    method = _get_transform_method(method, use_scipy)

    # Transform the image using the appropriate function
    rotated_image = _supported_rotation_methods[method](image, rmatrix, shift, order, missing, clip)

    return rotated_image


def _get_transform_method(method, use_scipy):
    # This is re-used in affine_transform and GenericMap.rotate
    if method not in _supported_rotation_methods:
        raise ValueError(f'Method {method} not in supported methods: '
                         f'{_supported_rotation_methods.keys()}')

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


def _handle_clipping_externally(rotation_function):
    """
    Decorator for a rotation function to handle clipping externally
    """
    @wraps(rotation_function)
    def wrapper(image, matrix, shift, order, missing, clip):
        # Get the rotated image without any clipping
        rotated_image = rotation_function(image, matrix, shift, order, missing, False)

        if clip and not np.all(np.isnan(rotated_image)):
            # Clip the image to the input range, with checking for whether `missing` was used
            lower = np.nanmin([np.nanmax([missing, np.nanmin(rotated_image)]), np.nanmin(image)])
            upper = np.nanmax([np.nanmin([missing, np.nanmax(rotated_image)]), np.nanmax(image)])
            rotated_image.clip(lower, upper, out=rotated_image)

        return rotated_image

    return wrapper


def _handle_image_nans_externally(rotation_function):
    """
    Decorator for a rotation function to handle any NaNs in the image externally to the function
    """
    @wraps(rotation_function)
    def wrapper(image, matrix, shift, order, missing, clip):
        # Check for any NaNs and set them to a value that is in the range of the input image
        isnan = np.isnan(image)
        has_nans = np.any(isnan)
        nanfree_image = np.nan_to_num(image, nan=np.nanmin(image)) if has_nans else image

        # Pass the NaN-free image on to the rotation function
        rotated_image = rotation_function(nanfree_image, matrix, shift, order, missing, clip)

        # Restore the NaNs
        if has_nans:
            # Use a convolution to find all pixels that are affected by NaNs
            # We want a kernel size that is an odd number that is at least order+1
            size = 2*int(np.ceil(order/2))+1
            expanded_nans = convolve2d(isnan.astype(int),
                                       np.ones((size, size)).astype(int),
                                       mode='same')
            rotated_nans = scipy.ndimage.affine_transform(expanded_nans.T, matrix,
                                                          offset=shift, order=1,
                                                          mode='nearest').T
            rotated_image[rotated_nans > 0] = np.nan

        return rotated_image

    return wrapper


# Define individual rotation implementations below and insert in the dictionary at the end.
#
# The arguments must be, in order:
#   image : `numpy.ndarray`
#       The image, which could be integers or floats, and can include NaNs
#   matrix : `numpy.ndarray` that is 2x2
#       The linear transformation matrix (e.g., rotation+scale+skew)
#   shift : 2-element `numpy.ndarray`
#       The translational shift to apply to the image in each axis
#   order : `int`
#       The numerical parameter that controls the degree of interpolation
#   missing : `float`
#       The value to use for outside the bounds of the original image, and may be NaN
#   clip : `bool`
#       Whether to clip the output image to the range of the input image
#
# There are two decorators that can be used:
#
# * `_handle_clipping_externally` is for rotation functions that do not implement clipping
# * `_handle_image_nans_externally` is for rotation functions that do not handle NaNs in the image
#
# Do not modify `image` directly; copy first if necessary.


@_handle_clipping_externally
@_handle_image_nans_externally
def _rotation_scipy(image, matrix, shift, order, missing, clip):
    """
    Rotation using SciPy
    """
    rotated_image = scipy.ndimage.affine_transform(image.T, matrix, offset=shift, order=order,
                                                   mode='constant', cval=missing).T

    return rotated_image


@_handle_clipping_externally
@_handle_image_nans_externally
def _rotation_skimage(image, matrix, shift, order, missing, clip):
    """
    Rotation using scikit-image

    Although :func:`skimage.transform.warp` implements clipping, its clipping
    behavior is inconsistent across interpolation orders, so we choose to handle
    clipping externally.

    :func:`skimage.transform.warp` does not currently support setting `missing` to
    NaN for interpolation orders 2, 4, or 5, so `missing` will be set to zero then.
    """
    import skimage.transform

    # Make the rotation matrix 3x3 to include translation of the image
    skmatrix = np.zeros((3, 3))
    skmatrix[:2, :2] = matrix
    skmatrix[2, 2] = 1.0
    skmatrix[:2, 2] = shift
    tform = skimage.transform.AffineTransform(skmatrix)

    if np.isnan(missing) and order in {2, 4, 5}:
        warn_user("Setting `missing` to NaN is not supported for this order of scikit-image "
                  "rotation, so using zero instead.")
        missing = 0

    if issubclass(image.dtype.type, numbers.Integral):
        warn_user("Integer input data has been cast to float64.")
        adjusted_image = image.astype(np.float64)
    else:
        adjusted_image = image.copy()

    # Scale image to range [0, 1]
    im_min = np.nanmin([missing, np.min(adjusted_image)])
    adjusted_image -= im_min
    adjusted_missing = missing - im_min
    im_max = np.nanmax([adjusted_missing, np.max(adjusted_image)])
    if im_max > 0:
        adjusted_image /= im_max
        adjusted_missing /= im_max

    rotated_image = skimage.transform.warp(adjusted_image, tform, order=order,
                                           mode='constant', cval=adjusted_missing, clip=clip)

    # Convert the image back to its original range
    if im_max > 0:
        rotated_image *= im_max
    rotated_image += im_min

    return rotated_image


_supported_rotation_methods = {'scipy': _rotation_scipy,
                               'skimage': _rotation_skimage}
