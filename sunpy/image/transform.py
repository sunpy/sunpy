"""
Functions for geometrical image transformation and warping.
"""
import numbers
import warnings
import types

import numpy as np
import scipy.ndimage.interpolation

from sunpy.util.exceptions import SunpyUserWarning

__all__ = ['affine_transform']

def affine_transform(image, rmatrix, order=3, scale=1.0, image_center=None,
                     recenter=False, missing=0.0, method='skimage'):
    """
    Rotates, shifts and scales an image.

    Will use one of `skimage.transform.warp`, `scipy.ndimage.interpolation.affine_transform`,
    or a passed-in custom method as selected by `method=`.
    If the appropriate library is not installed, will default
    to `scipy.ndimage.interpolation.affine_transform`.

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
    method : `string` or :func:, optional
        1. If `string`: `skimage` or `scipy`
        If `scipy`, uses :func:`scipy.ndimage.interpolation.affine_transform`.
        If `skimage`, uses :func:`skimage.transform.warp`.
        2. Elif :func:, uses user-defined function to perform affine transform.
        See `notes` for function requirements.
        Default: `skimage`; otherwise on ImportError, will use `scipy`.

    Returns
    -------
    `numpy.ndarray`:
        New rotated, scaled and translated image.

    Notes
    -----
    This algorithm uses an affine transformation as opposed to a polynomial
    geometrical transformation, which by default is `skimage.transform.warp`.
    One can specify using `scipy.ndimage.interpolation.affine_transform`  as an alternative
    affine transformation. The two transformations use different algorithms
    and thus do not give identical output.

    When using for `skimage.transform.warp` with order >= 4 or using
    `scipy.ndimage.interpolation.affine_transform` at all, "NaN" values will
    replaced with zero prior to rotation. No attempt is made to retain the NaN
    values.

    Input arrays with integer data are cast to float 64 and can be re-cast using
    `numpy.ndarray.astype` if desired.

    Although this function is analogous to the IDL's ``rot`` function, it does not
    use the same algorithm as the IDL ``rot`` function.
    IDL's ``rot`` calls the `POLY_2D <https://www.harrisgeospatial.com/docs/poly_2d.html>`__
    method to calculate the inverse mapping of original to target pixel
    coordinates. This is a polynomial geometrical transformation.
    Then optionally it uses a bicubic convolution interpolation
    algorithm to map the original to target pixel values.

    If a custom transform function is passed to `method`, it must have the function signature:
    `foo(image,rmatrix,order,scale,missing,image_center,recenter)`
    (identical to the parent affine_transform, without the `method` argument),
    and it must return the new rotated, scaled, translated image.

    """
    if isinstance(method, str):
        if method == 'skimage':
            method = _skimage_affine_transform
        elif method == 'scipy':
            method = _scipy_affine_transform
        else:
            raise ValueError("`method` {} not supported.".format(method))
    if not isinstance(method, types.FunctionType):
        raise ValueError("argument `method` must be a string or function")
    
    # do affine transform with selected method
    try:
        rotated_image = method(image=image, rmatrix=rmatrix,
                               order=order, scale=scale,
                               missing=missing, image_center=image_center,
                               recenter=recenter)
    except (ImportError, NameError):
        warnings.warn('library in affine transform `method` was not found. Using `scipy` method',
                      SunpyUserWarning)
        rotated_image = _scipy_affine_transform(image, rmatrix, order,
                                                scale, missing, image_center,
                                                recenter)
        
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
    rmatrix = rmatrix / scale
    shift = _calculate_shift(image, rmatrix, image_center, recenter)
    
    if np.any(np.isnan(image)):
        warnings.warn("Setting NaNs to 0 for SciPy rotation.", SunpyUserWarning)

    # Transform the image using the scipy affine transform
    return scipy.ndimage.interpolation.affine_transform(
        np.nan_to_num(image).T, rmatrix, offset=shift, order=order,
        mode='constant', cval=missing).T


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
