"""
Functions for geometrical image transformation and warping.
"""
import numbers
import warnings

import numpy as np
import scipy.ndimage.interpolation

from sunpy.util.exceptions import SunpyUserWarning

__all__ = ['affine_transform']

# Permitted argument values for affine_transform 'method='
allowed_methods = ['scipy', 'skimage', 'opencv']

# Flags for converting input order to appropriate cv2 interpolation flag
# As of Sept. 2020, OpenCV warpAffine does not support order 2,4,5
_CV_ORDER_FLAGS = {
    0: cv2.INTER_NEAREST,
    1: cv2.INTER_LINEAR,
    2: cv2.INTER_CUBIC,
    3: cv2.INTER_CUBIC,
    4: cv2.INTER_CUBIC
    5: cv2.INTER_CUBIC
}

def affine_transform(image, rmatrix, order=3, scale=1.0, image_center=None,
                     recenter=False, missing=0.0, use_scipy=False):
    """
    Rotates, shifts and scales an image.

    Will use `skimage.transform.warp` unless scikit-image can't be imported
    then it will use`scipy.ndimage.affine_transform`.

    Parameters
    ----------
    image : `numpy.ndarray`
        2D image to be rotated.
    rmatrix : `numpy.ndarray` that is 2x2
        Linear transformation rotation matrix.
    order : `int` 0-5, optional
        Interpolation order to be used, defaults to 3. When using scikit-image this parameter
        is passed into `skimage.transform.warp` (e.g., 3 corresponds to bi-cubic interpolation).
        When using scipy it is passed into
        `scipy.ndimage.affine_transform` where it controls the order of the spline.
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
    use_scipy : `bool`, optional
        Force use of `scipy.ndimage.affine_transform`.
        Will set all "NaNs" in image to zero before doing the transform.
        Defaults to `False`, unless scikit-image can't be imported.

    Returns
    -------
    `numpy.ndarray`:
        New rotated, scaled and translated image.

    Notes
    -----
    This algorithm uses an affine transformation as opposed to a polynomial
    geometrical transformation, which by default is `skimage.transform.warp`.
    One can specify using `scipy.ndimage.affine_transform` as
    an alternative affine transformation. The two transformations use different
    algorithms and thus do not give identical output.

    When using `skimage.transform.warp` with order >= 4 or using
    `scipy.ndimage.affine_transform` at all, "NaN" values will be replaced with
    zero prior to rotation. No attempt is made to retain the "NaN" values.

    Input arrays with integer data are cast to float 64 and can be re-cast using
    `numpy.ndarray.astype` if desired.

    In the case of `skimage.transform.warp`, the image is normalized to [0, 1]
    before passing it to the function. It is later rescaled back to the original range.

    Although this function is analogous to the IDL's ``rot`` function, it does not
    use the same algorithm as the IDL ``rot`` function.
    IDL's ``rot`` calls the `POLY_2D <https://www.harrisgeospatial.com/docs/poly_2d.html>`__
    method to calculate the inverse mapping of original to target pixel
    coordinates. This is a polynomial geometrical transformation.
    Then optionally it uses a bicubic convolution interpolation
    algorithm to map the original to target pixel values.
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
    if not use_scipy:
        try:
            import skimage.transform
        except ImportError:
            warnings.warn("scikit-image could not be imported. Image rotation will use scipy",
                          ImportWarning)
            use_scipy = True
    if use_scipy:
        if np.any(np.isnan(image)):
            warnings.warn("Setting NaNs to 0 for SciPy rotation.", SunpyUserWarning)
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

        if issubclass(image.dtype.type, numbers.Integral):
            warnings.warn("Integer input data has been cast to float64.",
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

        rotated_image = skimage.transform.warp(adjusted_image, tform, order=order,
                                               mode='constant', cval=adjusted_missing)

        # Convert the image back to its original range if it is valid
        if not is_nan_image:
            if im_max > 0:
                rotated_image *= im_max
            rotated_image += im_min

    return rotated_image

def affine_transform_new(image, rmatrix, order=3, scale=1.0, image_center=None,
                         recenter=False, missing=0.0, method='skimage'):
    """
    Rotates, shifts and scales an image.

    Will use one of `skimage.transform.warp`, `scipy.ndimage.affine_transform`, 
    or `cv2.warpAffine`, as selected by `method=`. If the appropriate library is not 
    installed, will default to `scipy.ndimage.affine_transform`.

    Parameters
    ----------
    image : `numpy.ndarray`
        2D image to be rotated.
    rmatrix : `numpy.ndarray` that is 2x2
        Linear transformation rotation matrix.
    order : `int` 0-5, optional
        Interpolation order to be used, defaults to 3. When using scikit-image this parameter
        is passed into `skimage.transform.warp` (e.g., 3 corresponds to bi-cubic interpolation).
        When using scipy it is passed into
        `scipy.ndimage.affine_transform` where it controls the order of the spline.
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
    method : `string`, optional
        `skimage`, `scipy`, or `opencv`: chosen library for the affine transform.
        Defaults to `skimage`; if selected library is not installed, reverts to `scipy`.

    Returns
    -------
    `numpy.ndarray`:
        New rotated, scaled and translated image.

    Notes
    -----
    This algorithm uses an affine transformation as opposed to a polynomial
    geometrical transformation, which by default is `skimage.transform.warp`.
    One can specify using `scipy.ndimage.affine_transform`  or 'cv2.warpAffine` as
    an alternative affine transformation. The three transformations use different
    algorithms and thus do not give identical output.

    When using for `skimage.transform.warp` with order >= 4 or using
    `scipy.ndimage.affine_transform` at all, "NaN" values will
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
    """
    if method.lower() not in allowed_methods:
        raise ValueError("Input method {} not found in allowed options: {}".format(
            method,allowed_methods))

    # test library import; default to scipy if needed
    if method == 'opencv':
        try:
            from cv2 import warpAffine
        except ImportError:
            warnings.warn("cv2.warpAffine could not be imported.",
                          "Image rotation will use scipy", ImportWarning)
            method = 'scipy'

    if method = 'skimage':
        try:
            import skimage.transform
        except ImportError:
            warnings.warn("scikit-image could not be imported. Image rotation will use scipy",
                          ImportWarning)
            method = 'scipy'

    
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
            
    # do affine transform with selected method
    if method == 'opencv':
        rotated_image = _opencv_affine_transform(image, rmatrix, order, scale, missing, shift)
    elif method == 'skimage':
        rotated_image = _skimage_affine_transform(image, rmatrix, order, scale, missing, shift)
    else:
        rotated_image = _scipy_affine_transform(image, rmatrix, order, scale, missing, shift)
    
    return rotated_image


def _opencv_affine_transform(image, rmatrix, order, scale, missing, shift):
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
    shift : `sequence`
        Offset into the image where the transform is applied. Should contain one
        value for each axis of `image`.
        
    NOTE: `cv2.warpAffine` uses a different coordinate system than its `skimage` 
    and `scipy` counterparts; this function is designed to take in arguments utilizing 
    the `skimage` and `scipy` reference frame and convert them to the `cv2` frame.
    """

    # needed to convert missing from np.dtype to native type
    # required for warpAffine input
    try:
        missing = missing.tolist()
    except AttributeError:
        pass

    # convert order to appropriate cv2 flag
    if order not in [0,1,3]:
        warnings.warn("Input order {} not supported. Setting order to 3 for openCV rotation.",
                      SunpyUserWarning)
    order = _CV_ORDER_FLAGS[max(0,min(order,5))]

    # get appropriate cv transform matrix
    # (convert default coordinate system to openCV coordinate system) 
    rmatrix = rmatrix * scale
    trans = np.eye(3,3)
    rot_scale = np.eye(3,3)
    trans[:2,2] = [-shift[0],-shift[1]]
    rot_scale[:2,:2] = rmatrix.T
    rmatrix = (rot_scale @ trans)[:2]

    if issubclass(image.dtype.type, numbers.Integral):
        warnings.warn("Integer input data has been cast to float64, "
                      "which is required for the opencv2-image transform.",
                      SunpyUserWarning)
        adjusted_image = image.astype(np.float64)
    else:
        adjusted_image = image.copy()
        
    h,w = adjusted_image.shape
    return cv2.warpAffine(adjusted_image, rmatrix, (w,h), flags=order,borderMode=cv2.BORDER_CONSTANT, borderValue=missing)

def _skimage_affine_transform(image, rmatrix, order, scale, missing, shift):
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
    shift : `sequence`
        Offset into the image where the transform is applied. Should contain one
        value for each axis of `image`. 
    """

    rmatrix = rmatrix / scale

    # Make the rotation matrix 3x3 to include translation of the image
    skmatrix = np.zeros((3, 3))
    skmatrix[:2, :2] = rmatrix
    skmatrix[2, 2] = 1.0
    skmatrix[:2, 2] = shift
    tform = skimage.transform.AffineTransform(skmatrix)

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
        
    return skimage.transform.warp(adjusted_image, tform, order=order,
                                           mode='constant', cval=missing,
                                           preserve_range=True)

def _scipy_affine_transform(image, rmatrix, order, scale, missing, shift):
    """
    Apply `scipy.ndimage.affine_transform` to `image`

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
    shift : `sequence`
        Offset into the image where the transform is applied. Should contain one
        value for each axis of `image`.
        
    """
    
    rmatrix = rmatrix / scale

    if np.any(np.isnan(image)):
        warnings.warn("Setting NaNs to 0 for SciPy rotation.", SunpyUserWarning)

    # Transform the image using the scipy affine transform
    return scipy.ndimage.interpolation.affine_transform(
        np.nan_to_num(image).T, rmatrix, offset=shift, order=order,
        mode='constant', cval=missing).T
