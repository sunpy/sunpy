"""Image rotation function(s)"""
from __future__ import absolute_import

import numpy as np
from scipy.ndimage import interpolation as sp
try:
    from skimage import transform as sk
    force_scipy = False
except:
    force_scipy = True

__all__ = ['affine_transform']

#TODO: Add functionality to specify interpolation method and missing value
def affine_transform(image, rmatrix, order=4, scale=1.0, rotation_center=None,
                     recenter=False, missing=0.0, scipy=False):
    """    
    Rotates and shifts an image using an affine transform. Intended to replace Map.rotate().
    Currently uses the old C extension function to rotate and shif the image, though this will be
    replaced with scikit-image's AffineTransform class and warp() function.

    Parameters
    ----------
    image: ndarray
        2D Image to be rotated.
    rmatrix: 2x2
        Linear transformation rotation matrix.
    order: int 0-5
        Interpolation order to be used, when using scikit-image this parameter
        is passed into :fun:`skimage.transform.warp`.
        When using scipy it is passed into 
        :fun:`scipy.ndimage.interpolation.affine_transform` where it controls 
        thwe order of the spline.
    scale: float
        A scale factor for the image. Default is no scaling.
    rotation_center: tuple
        The point in the image to rotate around (axis of rotation).
        Default: centre of the array.
    recenter: bool or array-like
        Move the axis of rotation to the centre of the array or recenter coords.
        Default: True, recentre to the centre of the array.
    missing: flot
        The value to replace any missing data after the transformation.
    scipy: bool
        Force use of :fun:`scipy.ndimage.interpolation.affine_transform`.
        This defaults to False unless sckit-image is not installed.

    Returns
    -------
    New rotated, scaled and translated image.
    """

    rmatrix = rmatrix / scale
    array_center = (np.array(image.shape)-1)/2.0
    # A rename to make things clearer. Also make sure it's an array
    # TODO: Deal with this properly and change it in the keywords
    if rotation_center is not None:
        image_center = np.asanyarray(rotation_center)
    else:
        image_center = array_center
    if recenter == True:
        rot_center = array_center
    elif recenter == False:
        rot_center = image_center

    displacement = np.dot(rmatrix, image_center)
    shift = rot_center - displacement

    if scipy or force_scipy:
        # This is the scipy call
        rotated_image = sp.affine_transform(image, rmatrix, offset=shift,
                                            order=order, mode='constant',
                                            cval=missing)
    else:
        skmatrix = np.zeros((3, 3))
        skmatrix[:2, :2] = rmatrix
        skmatrix[2, 2] = 1.0
        skmatrix[:2, 2] = shift
        im_max = image.max()
        tform = sk.AffineTransform(skmatrix)
        rotated_image = sk.warp(image, tform, order=order, mode='constant',
                                cval=missing) * im_max
    

    return rotated_image
