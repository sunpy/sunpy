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
def affine_transform(image, rmatrix=None, angle=None, scale=1.0,
                     rotation_center=None, recenter=True, scipy=False,
                     missing=0.0, interp_type='biquartic', interp_param=None):
    """    
    Rotates and shifts an image using an affine transform. Intended to replace Map.rotate().
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

    interps = ['nearest', 'bilinear', 'biquadratic', 'bicubic', 'biquartic',
               'biquintic', 'spline']
    if force_scipy or scipy:
        interp_type = 'spline'
    assert interp_type in interps
    if interp_param is None:
        if interp_type is 'spline':
            interp_param = 3
        elif interp_type is 'bicubic':
            interp_param = 0.5 # Should this be -0.5?
        else:
            interp_param = 0

    # A rename to make things clearer.
    # TODO: Deal with this properly and change it in the keywords
    image_center = rotation_center
    rmatrix = rmatrix / scale
    array_center = (np.array(image.shape)-1)/2.0
    if recenter == True:
        rot_center = array_center
    elif recenter == False:
        rot_center = image_center

    displacement = np.array([rmatrix[0,0]*image_center[1] + rmatrix[1,0]*image_center[0],
                        rmatrix[0,1]*image_center[1] + rmatrix[1,1]*image_center[0]])
    shift = rot_center - displacement

    if interp_type == 'spline':
        # This is the scipy call
        rotated_image = sp.affine_transform(image, rmatrix, offset=shift,
                                            order=interp_param, mode='constant',
                                            cval=missing)
    else:
        skmatrix = np.zeros((3, 3))
        skmatrix[:2, :2] = rmatrix
        skmatrix[2, 2] = 1.0
        skmatrix[:2, 2] = [shift[1], shift[0]]
        order = interps.index(interp_type)
        im_max = image.max()
        tform = sk.AffineTransform(skmatrix)
        rotated_image = sk.warp(image, tform, order=order, mode='constant',
                                cval=missing) * im_max
    

    return rotated_image
