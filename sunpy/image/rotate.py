"""Image rotation function(s)"""
from __future__ import absolute_import

import sunpy.image.Crotate as Crotate

__all__ = ['affine_transform']

#TODO: Add functionality to rescale image, specify whether or not to recentre, and specify rotation centre, interpolation method and missing value
def affine_transform(image, rmatrix=None, shift=None, angle=None):
    """Rotates and shifts an image using an affine transform. Intended to replace Map.rotate().
    Currently uses the old C extension function to rotate and shif the image, though this will be
    replaced with scikit-image's AffineTransform class and warp() function.

    Parameters
    ----------
    image: ndarray
        Image to be rotated.
    rmatrix: NxN
        Linear transformation rotation matrix. Either rmatrix or angle should be specified (not both).
    shift: tuple
        Amount by which to translate image (optional).
    angle: float
        Angle to rotate the image by (in radians). Either angle or rmatrix should be specified (not both).

    Returns
    -------
    Image rotated by angle or using rmatrix
    """

    rotated_image = Crotate.affine_transform(image, rmatrix, offset=shift,
                kernel=Crotate.BICUBIC, cubic=-0.5, mode='constant',
                cval=image.min())

    return rotated_image
