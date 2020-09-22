"""
Functions for geometrical image transformation and warping.
"""
import sys
import types
import numbers

import numpy as np

from sunpy.util.exceptions import warn_deprecated, warn_user

__all__ = ['add_rotation_function', 'affine_transform']


def affine_transform(image, rmatrix, order=3, scale=1.0, image_center=None,
                     recenter=False, missing=0.0, method='skimage', use_scipy=False):
    """
    Rotates, shifts and scales an image.

    Will use one of :func:`skimage.transform.warp`, :func:`scipy.ndimage.interpolation.affine_transform`,
    :func:`cv2.warpAffine`, or a passed-in custom method as selected by ``method=``.
    If the appropriate library is not installed, will raise ImportError.

    See `notes` for description of algorithm and definition of coordinate system.

    Parameters
    ----------
    image : `numpy.ndarray`
        2D image to be rotated.
    rmatrix : `numpy.ndarray` that is 2x2
        Linear transformation rotation matrix.
    order : `int` 0-5, optional
        Interpolation order to be used, defaults to 3. When using ``scikit-image`` this parameter
        is passed into :func:`skimage.transform.warp` (e.g., 3 corresponds to bi-cubic interpolation).
        When using ``scipy`` it is passed into :func:`scipy.ndimage.interpolation.affine_transform`
        where it controls the order of the spline.
        When using ``openCV``, it is converted into the appropriate flag; only order=0,1,3
        is supported; other values will raise an error
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
    method : `str` or function(), optional
        1. If a string:  one of ``'skimage'``, ``'scipy'``, or ``'cv2'``
        If ``'skimage'``, uses :func:`skimage.transform.warp`.
        If ``'scipy'``, uses :func:`scipy.ndimage.interpolation.affine_transform`.
        If ``'cv2'``, uses :func:`cv2.warpAffine`.
        2. If a function, uses the user-defined function to perform affine transform.
        See `notes` for function requirements.
        Default: ``'skimage'``: Will attempt to use :func:`skimage.transform.warp`;
        on ImportError, will use :func:`scipy.ndimage.interpolation.affine_transform`.
        (This behavior is identical to the now-deprecated ``use_scipy=False``)

    Returns
    -------
    `numpy.ndarray`:
        New rotated, scaled and translated image.

    Notes
    -----
    This algorithm uses an affine transformation as opposed to a polynomial
    geometrical transformation, which by default is :func:`skimage.transform.warp`.
    One can specify using :func:`scipy.ndimage.interpolation.affine_transform` or
    :func:`cv2.warpAffine` as alternative affine transformations. The transformations
    use different algorithms and thus do not give identical output.

    When using for :func:`skimage.transform.warp` with order >= 4 or using
    :func:`scipy.ndimage.interpolation.affine_transform` at all, NaN values will
    replaced with zero prior to rotation. No attempt is made to retain the NaN
    values.

    When using :func:`cv2.warpAffine`, only order=0,1,3 are supported. order=2,4,5 will
    be cast to order = 3 (bi-cubic interpolation).

    Input arrays with integer data are cast to float 64 and can be re-cast using
    `numpy.ndarray.astype` if desired.

    Although this function is analogous to the IDL's ``rot`` function, it does not
    use the same algorithm as the IDL ``rot`` function.
    IDL's ``rot`` calls the `POLY_2D <https://www.harrisgeospatial.com/docs/poly_2d.html>`__
    method to calculate the inverse mapping of original to target pixel
    coordinates. This is a polynomial geometrical transformation.
    Then optionally it uses a bicubic convolution interpolation
    algorithm to map the original to target pixel values.

    Method
    ------
    If a custom transform function is passed to ``method``, it must have the function signature:
    ``custom_transform(image, rmatrix, order, scale, missing, image_center, recenter)``
    (identical to the parent affine_transform, without the ``method`` argument),
    and it must return the new rotated, scaled, translated image.

    Each of the built-in methods are defined such that
    ``rmatrix`` = | cos(a) | -sin(a) | goes counterclockwise by an angle a
                  | sin(a) |  cos(a) |
    and the rotation axis is defined by ``image_center``.

    Since the built-in methods, by default, apply rotation about (0,0) (the upper left corner
    of the original image), a helper function ``_calculate_shift`` is provided
    to calculate an appropriate shift for the image such that the combined
    translation, rotation, and scaling yields the correct result.
    (see :func:`image.transform._calculate_shift` for details)
    ``shift`` = [a,b] such that the image is translated ``a`` pixels left, and ``b`` pixels up.
    """
    _allowed_methods = {'scipy': _scipy_affine_transform, 'skimage': _skimage_affine_transform,
                        'cv2': _opencv_affine_transform}

    if use_scipy:
        warn_deprecated("Argument `use_scipy` has been deprecated."
                        " Please set `method='scipy'` in the future.")
        method = 'scipy'

    if method == 'skimage':
        try:
            import skimage  # NoQA
        except ImportError:
            warn_user('`skimage` could not be imported. Using `scipy` instead')
            warn_deprecated('This fallback behavior will be deprecated. '
                            'Future versions will throw an ImportError and cease execution.')
            method = 'scipy'

    if isinstance(method, str):
        try:
            method = _allowed_methods[method]
        except KeyError:
            raise ValueError("`method` {} not supported.".format(method))

    if not isinstance(method, types.FunctionType):
        raise ValueError("argument `method` must be a string or function")
    # do affine transform with selected method
    try:
        rotated_image = method(image=image, rmatrix=rmatrix,
                               order=order, scale=scale,
                               missing=missing, image_center=image_center,
                               recenter=recenter)
    except ImportError:
        raise ImportError("Selected method not installed.")

    return rotated_image


def _skimage_affine_transform(image, rmatrix, order, scale, missing, image_center, recenter):
    """
    Apply :func:`skimage.transform.warp` to `image`

    Parameters
    ----------
    image: `numpy.ndarray`
        2D image to be rotated
    rmatrix : `numpy.ndarray` that is 2x2
        Linear transformation rotation matrix.
    order : `int` 0-5
        Interpolation order to be used
    scale : `float`
        A scale factor for the image
    missing : `float`
        The value to replace any missing data after the transformation.
    image_center : tuple, optional
        The point in the image to rotate around (axis of rotation).
        Defaults to the center of the array.
    recenter : `bool` or array-like, optional
        Move the axis of rotation to the center of the array or recenter coords.
        Defaults to `True` i.e., recenter to the center of the array.
    """
    from skimage.transform import AffineTransform, warp

    rmatrix = rmatrix / scale

    shift = _calculate_shift(image, rmatrix, image_center, recenter)

    # Make the rotation matrix 3x3 to include translation of the image
    skmatrix = np.zeros((3, 3))
    skmatrix[:2, :2] = rmatrix
    skmatrix[2, 2] = 1.0
    skmatrix[:2, 2] = shift
    tform = AffineTransform(skmatrix)

    if issubclass(image.dtype.type, numbers.Integral):
        warn_user("Integer input data has been cast to float64, "
                  "which is required for the scikit-image transform.")
        adjusted_image = image.astype(np.float64)
    else:
        adjusted_image = image.copy()

    if np.any(np.isnan(adjusted_image)) and order >= 4:
        warn_user("Setting NaNs to 0 for higher-order scikit-image rotation.")
        adjusted_image = np.nan_to_num(adjusted_image)

    # Scale image to range [0, 1] if it is valid (not made up entirely of NaNs)
    is_nan_image = np.all(np.isnan(adjusted_image))
    if is_nan_image:
        adjusted_missing = missing
    else:
        im_min = np.nanmin(adjusted_image)
        adjusted_image -= im_min
        im_max = np.nanmax(adjusted_image)
        if im_max > 0:
            adjusted_image /= im_max
            adjusted_missing = (missing - im_min) / im_max
        else:
            # The input array is all one value (aside from NaNs), so no scaling is needed
            adjusted_missing = missing - im_min

    rotated_image = warp(adjusted_image, tform, order=order,
                         mode='constant', cval=adjusted_missing)

    # Convert the image back to its original range if it is valid
    if not is_nan_image:
        if im_max > 0:
            rotated_image *= im_max
        rotated_image += im_min

    return rotated_image


def _scipy_affine_transform(image, rmatrix, order, scale, missing, image_center, recenter):
    """
    Apply :func:`scipy.ndimage.interpolation.affine_transform` to `image`

    Parameters
    ----------
    image: `numpy.ndarray`
        2D image to be rotated
    rmatrix : `numpy.ndarray` that is 2x2
        Linear transformation rotation matrix.
    order : `int` 0-5
        Interpolation order to be used
    scale : `float`
        A scale factor for the image
    missing : `float`
        The value to replace any missing data after the transformation.
    image_center : tuple, optional
        The point in the image to rotate around (axis of rotation).
        Defaults to the center of the array.
    recenter : `bool` or array-like, optional
        Move the axis of rotation to the center of the array or recenter coords.
        Defaults to `True` i.e., recenter to the center of the array.

    """
    import scipy.ndimage.interpolation

    rmatrix = rmatrix / scale
    shift = _calculate_shift(image, rmatrix, image_center, recenter)

    if np.any(np.isnan(image)):
        warn_user("Setting NaNs to 0 for SciPy rotation.")

    # Transform the image using the scipy affine transform
    return scipy.ndimage.interpolation.affine_transform(
        np.nan_to_num(image).T, rmatrix, offset=shift, order=order,
        mode='constant', cval=missing).T


def _opencv_affine_transform(image, rmatrix, order, scale, missing, image_center, recenter):
    """
    Apply :func:`cv2.warpAffine` to `image`

    Parameters
    ----------
    image: `numpy.ndarray`
        2D image to be rotated
    rmatrix : `numpy.ndarray` that is 2x2
        Linear transformation rotation matrix.
    order : `int` 0-5
        Interpolation order to be used
    scale : `float`
        A scale factor for the image
    missing : `float`
        The value to replace any missing data after the transformation.
    image_center : tuple, optional
        The point in the image to rotate around (axis of rotation).
        Defaults to the center of the array.
    recenter : `bool` or array-like, optional
        Move the axis of rotation to the center of the array or recenter coords.
        Defaults to `True` i.e., recenter to the center of the array.

    """
    import cv2

    # Flags for converting input order from `integer` to the appropriate interpolation flag
    # As of Oct. 2021, OpenCV warpAffine does not support order 2,4,5
    _CV_ORDER_FLAGS = {
        0: cv2.INTER_NEAREST,
        1: cv2.INTER_LINEAR,
        3: cv2.INTER_CUBIC,
    }

    try:
        order = _CV_ORDER_FLAGS[order]
    except KeyError:
        raise ValueError("Input order={} not supported in openCV. ".format(order),
                         "Please use order=0, 1, or 3.")

    # needed to convert ``missing`` from potentially a np.dtype
    # to the native `int` type required for :func:`cv2.warpAffine`
    try:
        missing = missing.tolist()
    except AttributeError:
        pass

    # OpenCV applies the shift+rotation operations in a different order; we need to calculate
    # translation using ``rmatrix/scale``, but scale+rotation with ``rmatrix*scale``
    # in order to match what skimage/scipy do
    shift = _calculate_shift(image, rmatrix / scale, image_center, recenter)

    # get appropriate cv transform matrix
    # (with a slight amount of voodoo to adjust for different coordinate systems)
    rmatrix = rmatrix * scale

    trans = np.eye(3, 3)
    rot_scale = np.eye(3, 3)
    trans[:2, 2] = [-shift[0], -shift[1]]
    rot_scale[:2, :2] = rmatrix.T
    rmatrix = (rot_scale @ trans)[:2]

    if issubclass(image.dtype.type, numbers.Integral):
        warn_user("Integer input data has been cast to float64, "
                  "for the openCV transform.")
        adjusted_image = image.astype(np.float64)
    else:
        adjusted_image = image.copy()

    h, w = adjusted_image.shape

    # equivalent to skimage.transform.warp(adjusted_image, tform, order=order,
    #                                     mode='constant', cval=missing)
    return cv2.warpAffine(adjusted_image, rmatrix, (w, h), flags=order,
                          borderMode=cv2.BORDER_CONSTANT, borderValue=missing)


def _calculate_shift(image, rmatrix, image_center=None, recenter=False):
    """
    Returns the offset into the array where the transform is applied.

    Parameters
    ----------
    image: `numpy.ndarray`
        2D image to be rotated
    rmatrix : `numpy.ndarray` that is 2x2
        Linear transformation rotation matrix.
    image_center : tuple, optional
        The point in the image to rotate around (axis of rotation).
        Defaults to the center of the array.
    recenter : `bool` or array-like, optional
        Move the axis of rotation to the center of the array or recenter coords.
        Defaults to `True` i.e., recenter to the center of the array.

    Returns
    -------
    tuple :
       Offset into array where transform is applied, based on rotation and image center.

    Notes
    -----
    Assumptions:
    1. Rotation, as defined in rmatrix, is counterclockwise.
    i.e. |cos(angle) -sin(angle)|
         |sin(angle)  cos(angle)|
    2. Rotation matrix has already been scaled by a scale factor, if applicable
    (NOTE: the appropriate scaling is usually rmatrix/scale)
    3. Coordinate system of offset array ``[a,b]`` is such that image is shifted
    ``a`` pixels to the left and ``b`` pixels up.
    """

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
    * An input image with byte ordering that does not match the native byte order of
      the system (e.g., big-endian values on a little-endian system) will be
      copied and byte-swapped prior to rotation.
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

    # Swap the byte order if it is non-native (e.g., big-endian on a little-endian system)
    if adjusted_image.dtype.byteorder == ('>' if sys.byteorder == 'little' else '<'):
        adjusted_image = adjusted_image.byteswap().newbyteorder()

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
    return image_center - displacement
