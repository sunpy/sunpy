"""
Helper functions for image manipulation
"""
from __future__ import absolute_import, division

import numpy as np
from skimage.util import img_as_float

__all__ = ['to_norm', 'un_norm']


def to_norm(arr):
    """
    Helper function to normalise/scale an array.  This is needed for example
    for scikit-image which uses floats between 0 and 1.

    Parameters
    ----------
    arr : `~numpy.ndarray`
        Array to normalise.

    Returns
    -------
    arr : `~numpy.ndarray`
        Array with values between 0 (min) and 1 (max)

    Examples
    --------
    >>> import numpy as np
    >>> from sunpy.image.util import to_norm
    >>> out = to_norm(np.array([-1, 0, 1]))
    >>> out
    array([0. , 0.5, 1. ])
    """
    arr = np.array(arr, dtype='double')
    arr = img_as_float(arr, force_copy=True)
    if arr.min() < 0:
        arr += np.abs(arr.min())
    arr /= arr.max()
    return arr


def un_norm(arr, original):
    """
    Helper function to un-normalise (or re-scale) an array based in
    the values of the original array.

    Parameters
    ----------
    arr : `~numpy.ndarray`
        Array of floats to un-normalise with values in [0,1]
    original : `~numpy.ndarray`
        Original array with the min and max values

    Returns
    -------
    arr : `~numpy.ndarray`
        Array with values between `original.min()` and `original.max()` . Note
        that the type of the original image is not guaranteed to be reproduced.

    Examples
    --------
    >>> import numpy as np
    >>> from sunpy.image.util import un_norm
    >>> original = np.array([-1, 0, 1])
    >>> normalised = np.array([0., 0.5, 1.])
    >>> out = un_norm(normalised, original)
    >>> out
    array([-1.,  0.,  1.])
    """
    level = 0 if original.min() > 0 else np.abs(original.min())
    arr *= original.max() + level
    arr -= level
    return arr
