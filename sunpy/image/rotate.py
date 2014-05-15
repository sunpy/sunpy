"""Image rotation function(s)"""
from __future__ import absolute_import

import numpy as np
import sunpy.image.Crotate as Crotate
from skimage import transform as tf

__all__ = ['affine_transform']

#TODO: Add functionality to specify interpolation method and missing value
def affine_transform(image, rmatrix=None, angle=None, scale=1.0, rotation_center=None,
                     recenter=True, rotate_func='skimage'):
    """Rotates and shifts an image using an affine transform. Intended to replace Map.rotate().
    Currently uses the old C extension function to rotate and shif the image, though this will be
    replaced with scikit-image's AffineTransform class and warp() function.

    Parameters
    ----------
    image: ndarray
        Image to be rotated.
    rmatrix: NxN
        Linear transformation rotation matrix. Either rmatrix or angle should be specified (not both).
    angle: float
        Angle to rotate the image by (in radians). Either angle or rmatrix should be specified (not both).
    scale: float
        A scale factor for the image. Default is no scaling.
    rotation_center: tuple
        The point in the image to rotate around (axis of rotation).
        Default: centre of the array.
    recenter: bool or array-like
        Move the axis of rotation to the centre of the array or recenter coords.
        Default: True, recentre to the centre of the array.

    Returns
    -------
    New rotated, scaled and translated image.
    """

    assert rotate_func in ['skimage', 'Crotate']

    rmatrix = rmatrix / scale
    center = (np.array(image.shape)-1)/2.0
    if recenter == True:
        recenter = center
    elif recenter == False:
        recenter = rotation_center

    if rotation_center:
        shift = np.array(rotation_center) - np.array(recenter)
    else:
        shift = np.array([0.0, 0.0])

    if rotate_func == 'skimage':
        skmatrix = np.zeros((3, 3))
        skmatrix[:2, :2] = rmatrix
        skmatrix[2, 2] = 1.0
        skmatrix[:2, 2] = [shift[1], shift[0]]
        im_max = image.max()
        tform = tf.AffineTransform(skmatrix)
        rotated_image = tf.warp(image/im_max, tform, order=3,
                     mode='constant', cval=image.min()) * im_max
    elif rotate_func == 'Crotate':
        rotated_image = Crotate.affine_transform(image, rmatrix, offset=shift,
                    kernel=Crotate.BICUBIC, cubic=-0.5, mode='constant',
                    cval=image.min())

    return rotated_image
