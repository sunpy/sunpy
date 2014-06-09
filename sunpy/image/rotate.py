"""Image rotation function(s)"""
from __future__ import absolute_import

import numpy as np
import scipy.ndimage.interpolation
try:
    import skimage.transform
    scikit_image_not_found = False
except ImportError:  # pragma: no cover
    scikit_image_not_found = True  # pragma: no cover

__all__ = ['affine_transform']

def affine_transform(image, rmatrix, order=4, scale=1.0, image_center=None,
                     recenter=False, missing=0.0, use_scipy=False):
    """    
    Rotates, shifts and scales an image using :func:`skimage.transform.warp`, or
    :func:`scipy.ndimage.interpolation.affine_transform` if specified. Falls back to
    the scipy function if scikit-image can't be imported.

    Parameters
    ----------
    image: ndarray
        2D Image to be rotated.
    rmatrix: 2x2
        Linear transformation rotation matrix.
    order: int 0-5
        Interpolation order to be used. When using scikit-image this parameter
        is passed into :func:`skimage.transform.warp`.
        When using scipy it is passed into 
        :func:`scipy.ndimage.interpolation.affine_transform` where it controls 
        the order of the spline.
    scale: float
        A scale factor for the image. Default is no scaling.
    rotation_center: tuple
        The point in the image to rotate around (axis of rotation).
        Default: center of the array.
    recenter: bool or array-like
        Move the axis of rotation to the center of the array or recenter coords.
        Default: True, recenter to the center of the array.
    missing: float
        The value to replace any missing data after the transformation.
    scipy: bool
        Force use of :func:`scipy.ndimage.interpolation.affine_transform`.
        Default: False unless sckit-image is not installed.

    Returns
    -------
    New rotated, scaled and translated image.
    """

    rmatrix = rmatrix / scale
    array_center = (np.array(image.shape)-1)/2.0
    
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

    displacement = np.dot(rmatrix, image_center)
    shift = rot_center - displacement

    if use_scipy or scikit_image_not_found:
        # Transform the image using the scipy affine transform
        rotated_image = scipy.ndimage.interpolation.affine_transform(
                image, rmatrix, offset=shift, order=order, mode='constant',
                cval=missing)
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
        rotated_image = skimage.transform.warp(image/abs(image).max(), tform,
                                               order=order, mode='constant',
                                               cval=missing) * abs(image).max()

    return rotated_image
