"""
Functions for geometrical image transformation and warping.
"""
import warnings

import numpy as np
import scipy.ndimage.interpolation

from sunpy.util.exceptions import SunpyUserWarning

try:
    import skimage.transform
    scikit_image_not_found = False
except ImportError:
    warnings.warn("scikit-image could not be imported. Image rotation will use scipy",
                  ImportWarning)
    scikit_image_not_found = True


__all__ = ['affine_transform']


def affine_transform(image, rmatrix, order=3, scale=1.0, image_center=None,
                     recenter=False, missing=0.0, use_scipy=False):
    """
    Rotates, shifts and scales an image.

    Will use `skimage.transform.warp` unless scikit-image can't be imported
    then it will use`scipy.ndimage.interpolation.affine_transform`.

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
        `scipy.ndimage.interpolation.affine_transform` where it controls the order of the spline.
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
        Force use of `scipy.ndimage.interpolation.affine_transform`.
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
    One can specify using `scipy.ndimage.interpolation.affine_transform` as
    an alternative affine transformation. The two transformations use different
    algorithms and thus do not give identical output.

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

        # Transform the image using the skimage function
        if not np.issubdtype(image.dtype, np.float64):
            warnings.warn("Input data has been cast to float64.", SunpyUserWarning)
            adjusted_image = image.astype(np.float64)
        else:
            adjusted_image = image.copy()
        if np.any(np.isnan(adjusted_image)) and order >= 4:
            warnings.warn("Setting NaNs to 0 for higher-order scikit-image rotation.", SunpyUserWarning)
            adjusted_image = np.nan_to_num(adjusted_image)

        rotated_image = skimage.transform.warp(adjusted_image, tform, order=order,
                                               mode='constant', cval=missing)

    return rotated_image
