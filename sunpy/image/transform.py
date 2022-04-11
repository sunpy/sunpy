"""
Functions for geometrical image transformation and warping.
"""
import numbers
from functools import wraps

import numpy as np
import scipy.ndimage
from scipy.signal import convolve2d

from sunpy.util.exceptions import warn_deprecated, warn_user

__all__ = ['add_rotation_function', 'affine_transform']


def affine_transform(image, rmatrix, order=3, scale=1.0, image_center=None,
                     recenter=False, missing=np.nan, use_scipy=None, *, method='scipy', clip=True):
    """
    Rotates, shifts and scales an image.

    This function supports NaN values in the input image and supports using NaN for
    pixels in the output image that are beyond the extent of the input image.

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
        Defaults to `False`.
    missing : `float`, optional
        The value to use for pixels in the output image that are beyond the extent
        of the input image.
        Defaults to `numpy.nan`.
    method : {{{rotation_function_names}}}, optional
        Rotation function to use.
        Defaults to ``'scipy'``.
    clip : `bool`, optional
        If `True`, clips the pixel values of the output image to the range of the
        input image (including the value of ``missing``, if used).
        Defaults to `True`.

    Returns
    -------
    `numpy.ndarray`:
        New rotated, scaled and translated image.

    Notes
    -----
    For each NaN pixel in the input image, one or more pixels in the output image
    will be set to NaN, with the number of pixels affected depending on the
    interpolation order.

    See :func:`~sunpy.image.transform.add_rotation_function` for how to add a
    different rotation function.
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
    rotated_image = _rotation_function_registry[method](image, rmatrix, shift, order, missing, clip)

    return rotated_image


def _get_transform_method(method, use_scipy):
    # This is re-used in affine_transform and GenericMap.rotate
    if method not in _rotation_function_registry:
        raise ValueError(f'Method {method} not in supported methods: '
                         f'{_rotation_function_registry.keys()}')

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


def add_rotation_function(name, handles_clipping, handles_image_nans, handles_nan_missing):
    """
    Decorator to add a rotation function to the registry of selectable
    implementations.

    Each registered rotation function becomes a selectable option for
    :func:`sunpy.image.transform.affine_transform` and
    :meth:`sunpy.map.GenericMap.rotate`.  Those two routines are required to handle
    clipping the output image, NaNs in the input image, and NaN as the value to use
    for pixels in the output image that are beyond the extent of the input image.
    If the supplied rotation function cannot provide one or more of these capabilities,
    the decorator is able to provide them instead.

    The decorator accepts the parameters listed under ``Parameters``.  The decorated
    rotation function must accept the parameters listed under ``Other Parameters``
    in that order and return the rotated image.

    Parameters
    ----------
    name : `str`
        The name that will be used to select the rotation function
    handles_clipping : `bool`
        Specifies whether the rotation function will internally perform clipping.
        If ``False``, the rotation function will always receive ``False`` for the
        ``clip`` input parameter.
    handles_image_nans : `bool`
        Specifies whether the rotation function will internally handle NaNs in the
        input image.  If ``False``, the rotation function is guaranteed to be
        provided an image without any NaNs.
    handles_nan_missing : `bool`
        Specifies whether the rotation function will internally handle NaN as the
        ``missing`` value.  If ``False``, the rotation function will never receive
        NaN, but instead receive a value in the input range of the image.

    Other Parameters
    ----------------
    image : `numpy.ndarray`
        The image, which could be integers or floats
    matrix : `numpy.ndarray` that is 2x2
        The linear transformation matrix (e.g., rotation+scale+skew)
    shift : 2-element `numpy.ndarray`
        The translational shift to apply to the image in each axis
    order : `int`
        The numerical parameter that controls the degree of interpolation
    missing : `float`
        The value to use for outside the bounds of the original image
    clip : `bool`
        Whether to clip the output image to the range of the input image

    Notes
    -----
    The docstring of the rotation function should be a bulleted list of notes
    specific to the rotation function.  It will be appended to ``Notes`` section of
    the docstring for :func:`~sunpy.image.transform.affine_transform`.

    The rotation function is supplied the input image directly, so the function
    should not modify the image in place.

    Setting any of the ``handles_*`` parameters to ``False`` means that computation
    will be performed to modify the image returned by the rotation function before
    it is returned to :func:`~sunpy.image.transform.affine_transform`.

    If the decorator is handling image NaNs on behalf of the rotation function
    (i.e., ``handles_image_nans=False``), pixels in the output image will be set to
    NaN if they are within a number of pixels equal to half of the ``order``
    parameter.  This step requires an additional image convolution, which might be
    avoidable if the rotation function were able to internally handle image NaNs.
    """
    def decorator(rotation_function):
        @wraps(rotation_function)
        def wrapper(image, matrix, shift, order, missing, clip):
            clip_to_use = clip if handles_clipping else False

            # If missing cannot be used directly, use a value in the input range of the image
            needs_missing_handling = not handles_nan_missing and np.isnan(missing)
            missing_to_use = missing if not needs_missing_handling else np.nanmin(image)

            if not handles_image_nans:
                isnan = np.isnan(image)
                needs_nan_handling = np.any(isnan)
            else:
                needs_nan_handling = False

            # If needed, set any image NaNs to a value that is in the range of the input image
            image_to_use = np.nan_to_num(image, nan=np.nanmin(image)) if needs_nan_handling else image

            rotated_image = rotation_function(image_to_use, matrix, shift, order,
                                              missing_to_use, clip_to_use)

            # If needed, restore the NaNs
            if needs_nan_handling:
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

            if needs_missing_handling:
                # Rotate a constant image to determine where to apply the NaNs for `missing`
                constant = scipy.ndimage.affine_transform(np.ones_like(image).T.astype(int), matrix,
                                                          offset=shift, order=0, mode='constant').T
                rotated_image[constant < 1] = np.nan

            if not handles_clipping and clip and not np.all(np.isnan(rotated_image)):
                # Clip the image to the input range
                if np.isnan(missing):
                    # If `missing` is NaN, clipping to the input range is straightforward
                    rotated_image.clip(np.nanmin(image), np.nanmax(image), out=rotated_image)
                else:
                    # Otherwise, check if `missing` should be considered part of the input range
                    lower = np.nanmin([np.max([missing, np.nanmin(rotated_image)]), np.nanmin(image)])
                    upper = np.nanmax([np.min([missing, np.nanmax(rotated_image)]), np.nanmax(image)])
                    rotated_image.clip(lower, upper, out=rotated_image)

            return rotated_image

        _rotation_function_registry[name] = wrapper

        # Add the docstring of the rotation function to the docstring of affine_transform
        affine_transform.__doc__ += (f"\n    **Specific notes for the '{name}' rotation method:**"
                                     f"\n{wrapper.__doc__}")

        return wrapper
    return decorator


_rotation_function_registry = {}


@add_rotation_function("scipy",
                       handles_clipping=False, handles_image_nans=False, handles_nan_missing=True)
def _rotation_scipy(image, matrix, shift, order, missing, clip):
    """
    * Rotates using :func:`scipy.ndimage.affine_transform`
    * The ``mode`` parameter for :func:`~scipy.ndimage.affine_transform` is fixed to
      be ``'constant'``
    """
    rotated_image = scipy.ndimage.affine_transform(image.T, matrix, offset=shift, order=order,
                                                   mode='constant', cval=missing).T

    return rotated_image


@add_rotation_function("skimage",
                       handles_clipping=False, handles_image_nans=False, handles_nan_missing=False)
def _rotation_skimage(image, matrix, shift, order, missing, clip):
    """
    * Rotates using :func:`skimage.transform.warp`
    * The implementation for higher orders of interpolation means that the pixels
      in the output image that are beyond the extent of the input image may not have
      exactly the value of the ``missing`` parameter.
    * An input image with integer data is cast to floats prior to passing to
      :func:`~skimage.transform.warp`.  The output image can be re-cast using
      :meth:`numpy.ndarray.astype` if desired.
    * Does not let :func:`~skimage.transform.warp` handle clipping due to
      inconsistent handling across interpolation orders
    * Does not let :func:`~skimage.transform.warp` handle image NaNs because they
      are not handled properly for some interpolation orders
    * Does not pass NaN as ``missing`` to :func:`~skimage.transform.warp` due to
      inconsistent handling across interpolation orders
    """
    import skimage.transform

    # Make the rotation matrix 3x3 to include translation of the image
    skmatrix = np.zeros((3, 3))
    skmatrix[:2, :2] = matrix
    skmatrix[2, 2] = 1.0
    skmatrix[:2, 2] = shift
    tform = skimage.transform.AffineTransform(skmatrix)

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

    # Be aware that even though mode is set to 'constant', when skimage 0.19 calls scipy,
    # it specifies the scipy mode to be 'grid-constant' rather than 'constant'
    rotated_image = skimage.transform.warp(adjusted_image, tform, order=order,
                                           mode='constant', cval=adjusted_missing, clip=clip)

    # Convert the image back to its original range
    if im_max > 0:
        rotated_image *= im_max
    rotated_image += im_min

    return rotated_image


# Generate the string with allowable rotation-function names for use in docstrings
_rotation_function_names = ", ".join([f"``'{name}'``" for name in _rotation_function_registry])
# Insert into the docstring for affine_transform.  We cannot use the add_common_docstring decorator
# due to what would be a circular loop in definitions.
affine_transform.__doc__ = affine_transform.__doc__.format(rotation_function_names=_rotation_function_names)
