"""
Functions for geometrical image transformation and warping.
"""
import numbers
import warnings
import types

import numpy as np

from sunpy.util.exceptions import SunpyUserWarning, SunpyDeprecationWarning

__all__ = ['affine_transform']

def affine_transform(image, rmatrix, order=3, scale=1.0, image_center=None,
                     recenter=False, missing=0.0, method='skimage', use_scipy=False):
    """
    Rotates, shifts and scales an image.

    Will use one of `skimage.transform.warp`, `scipy.ndimage.interpolation.affine_transform`,
    `cv2.warpAffine`, or a passed-in custom method as selected by `method=`.
    If the appropriate library is not installed, will raise ImportError.

    See `notes` for description of algorithm and definition of coordinate system.

    Parameters
    ----------
    image : `numpy.ndarray`
        2D image to be rotated.
    rmatrix : `numpy.ndarray` that is 2x2
        Linear transformation rotation matrix.
    order : `int` 0-5, optional
        Interpolation order to be used, defaults to 3. When using scikit-image this parameter
        is passed into `skimage.transform.warp` (e.g., 3 corresponds to bi-cubic interpolation).
        When using scipy it is passed into`scipy.ndimage.interpolation.affine_transform`
        where it controls the order of the spline.
        When using openCV, it is converted into the appropriate flag; only order=0,1,3
        is supported; order=2,4,5 will be converted to bi-cubic (order = 3) interpolation.
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
    method : `string` or function(), optional
        1. If `string`: `skimage`, `scipy`, or `cv2'
        If `skimage`, uses :func:`skimage.transform.warp`.
        If `scipy`, uses :func:`scipy.ndimage.interpolation.affine_transform`.
        If `cv2`, uses :func:`cv2.warpAffine`.
        2. Elif function, uses user-defined function to perform affine transform.
        See `notes` for function requirements.
        Default: `skimage`: Will attempt to use :func:`skimage.transform.warp`;
        on ImportError, will use :func:`scipy.ndimage.interpolation.affine_transform`.
        (This behavior is identical to the now-deprecated `use_scipy=False`)

    Returns
    -------
    `numpy.ndarray`:
        New rotated, scaled and translated image.

    Notes
    -----
    This algorithm uses an affine transformation as opposed to a polynomial
    geometrical transformation, which by default is `skimage.transform.warp`.
    One can specify using `scipy.ndimage.interpolation.affine_transform` or
    `cv2.warpAffine` as alternative affine transformations. The transformations
    use different algorithms and thus do not give identical output.

    When using for `skimage.transform.warp` with order >= 4 or using
    `scipy.ndimage.interpolation.affine_transform` at all, "NaN" values will
    replaced with zero prior to rotation. No attempt is made to retain the NaN
    values.

    When using `cv2.warpAffine`, only order=0,1,3 are supported. order=2,4,5 will
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
    If a custom transform function is passed to `method`, it must have the function signature:
    `custom_transform(image,rmatrix,order,scale,missing,image_center,recenter)`
    (identical to the parent affine_transform, without the `method` argument),
    and it must return the new rotated, scaled, translated image.

    Each of the built-in methods are defined such that
    `rmatrix` = | cos(a) | -sin(a) | goes counterclockwise by an angle a
                | sin(a) |  cos(a) |
    and the rotation axis is defined by `image_center`.

    Since the built-in methods, by default, apply rotation about (0,0) (the upper left corner
    of the original image), a helper function `_calculate_shift` is provided
    to calculate an appropriate shift for the image such that the combined
    translation, rotation, and scaling yields the correct result.
    (see :func:`image.transform._calculate_shift` for details)
    `shift` = [a,b] such that the image is translated `a` pixels left, and `b` pixels up.
    """
    _allowed_methods = {'scipy': _scipy_affine_transform, 'skimage': _skimage_affine_transform,
                        'cv2': _opencv_affine_transform}

    if use_scipy:
        warnings.warn("Argument `use_scipy` has been deprecated."
                      " Please set `method='scipy'` in the future.",
                      SunpyDeprecationWarning)
        method = 'scipy'

    if method == 'skimage':
        try:
            import skimage
        except ImportError:
            warnings.warn('`skimage` could not be imported. Using `scipy` instead',
                          SunpyUserWarning)
            warnings.warn('This fallback behavior will be deprecated. '
                          'Future versions will throw an ImportError and cease execution.',
                          SunpyDeprecationWarning)
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
    Apply `skimage.transform.warp` to `image`

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
        warnings.warn("Integer input data has been cast to float64, "
                      "which is required for the scikit-image transform.",
                      SunpyUserWarning)
        adjusted_image = image.astype(np.float64)
    else:
        adjusted_image = image.copy()

    if np.any(np.isnan(adjusted_image)) and order >= 4:
        warnings.warn("Setting NaNs to 0 for higher-order scikit-image rotation.",
                      SunpyUserWarning)
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
    Apply `scipy.ndimage.interpolation.affine_transform` to `image`

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
        warnings.warn("Setting NaNs to 0 for SciPy rotation.", SunpyUserWarning)

    # Transform the image using the scipy affine transform
    return scipy.ndimage.interpolation.affine_transform(
        np.nan_to_num(image).T, rmatrix, offset=shift, order=order,
        mode='constant', cval=missing).T


def _opencv_affine_transform(image, rmatrix, order, scale, missing, image_center, recenter):
    """
    Apply cv2.warpAffine to `image`

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
    # As of Sept. 2020, OpenCV warpAffine does not support order 2,4,5
    _CV_ORDER_FLAGS = {
        0: cv2.INTER_NEAREST,
        1: cv2.INTER_LINEAR,
        3: cv2.INTER_CUBIC,
    }

    # convert order to appropriate cv2 flag
    try:
        order = _CV_ORDER_FLAGS[order]
    except KeyError:
        warnings.warn("input order={} not supported in openCV."
                      " order has been cast to 3".format(order),
                      SunpyUserWarning)
        order = _CV_ORDER_FLAGS[3]

    # needed to convert `missing` from potentially a np.dtype
    # to the native `int` type required for cv2.warpAffine
    try:
        missing = missing.tolist()
    except AttributeError:
        pass

    # OpenCV applies the shift+rotation operations in a different order; we need to calculate
    # translation using `rmatrix/scale`, but scale+rotation with `rmatrix*scale`
    # in order to match what skimage/scipy do
    shift = _calculate_shift(image, rmatrix/scale, image_center, recenter)

    # get appropriate cv transform matrix
    # (with a slight amount of voodoo to adjust for different coordinate systems)
    rmatrix = rmatrix*scale

    trans = np.eye(3, 3)
    rot_scale = np.eye(3, 3)
    trans[:2, 2] = [-shift[0], -shift[1]]
    rot_scale[:2, :2] = rmatrix.T
    rmatrix = (rot_scale @ trans)[:2]

    if issubclass(image.dtype.type, numbers.Integral):
        warnings.warn("Integer input data has been cast to float64, "
                      "for the openCV transform.",
                      SunpyUserWarning)
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
    3. Coordinate system of offset array [a,b] is such that image is shifted
    `a` pixels to the left and `b` pixels up.
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
    return image_center - displacement
