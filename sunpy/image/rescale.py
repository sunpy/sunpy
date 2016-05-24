"""Image resampling methods"""
from __future__ import absolute_import, division, print_function

import numpy as np
import scipy.interpolate
import scipy.ndimage
from sunpy.extern.six.moves import range

__all__ = ['resample', 'reshape_image_to_4d_superpixel']

def resample(orig, dimensions, method='linear', center=False, minusone=False):
    """Returns a new `numpy.ndarray` that has been resampled up or down.

    Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL's congrid routine (which apparently originally came from a
    VAX/VMS routine of the same name.)

    Parameters
    ----------
    orig : `~numpy.ndarray`
        Original inout array.
    dimensions : tuple
        Dimensions that new `~numpy.ndarray` should have.
    method : {'neighbor' | 'nearest' | 'linear' | 'spline'}
        Method to use for resampling interpolation.
            * neighbor - Closest value from original data
            * nearest and linear - Uses n x 1-D interpolations calculated by
              `scipy.interpolate.interp1d`.
            * spline - Uses ndimage.map_coordinates
    center : bool
        If True, interpolation points are at the centers of the bins,
        otherwise points are at the front edge of the bin.
    minusone : bool
        For orig.shape = (i,j) & new dimensions = (x,y), if set to False
        orig is resampled by factors of (i/x) * (j/y), otherwise orig
        is resampled by(i-1)/(x-1) * (j-1)/(y-1)
        This prevents extrapolation one element beyond bounds of input
        array.

    Returns
    -------
    out : `~numpy.ndarray`
        A new `~numpy.ndarray` which has been resampled to the desired
        dimensions.

    References
    ----------
    | http://www.scipy.org/Cookbook/Rebinning (Original source, 2011/11/19)
    """

    # Verify that number dimensions requested matches original shape
    if len(dimensions) != orig.ndim:
        raise UnequalNumDimensions("Number of dimensions must remain the same "
                                   "when calling resample.")

    #@note: will this be okay for integer (e.g. JPEG 2000) data?
    if orig.dtype not in [np.float64, np.float32]:
        orig = orig.astype(np.float64)

    dimensions = np.asarray(dimensions, dtype=np.float64)
    m1 = np.array(minusone, dtype=np.int64) # array(0) or array(1)
    offset = np.float64(center * 0.5)       # float64(0.) or float64(0.5)

    # Resample data
    if method == 'neighbor':
        data = _resample_neighbor(orig, dimensions, offset, m1)
    elif method in ['nearest','linear']:
        data = _resample_nearest_linear(orig, dimensions, method,
                                             offset, m1)
    elif method == 'spline':
        data = _resample_spline(orig, dimensions, offset, m1)
    else:
        raise UnrecognizedInterpolationMethod("Unrecognized interpolation "
                                              "method requested.")

    return data


def _resample_nearest_linear(orig, dimensions, method, offset, m1):
    """Resample Map using either linear or nearest interpolation."""

    dimlist = []

    # calculate new dims
    for i in range(orig.ndim):
        base = np.arange(dimensions[i])
        dimlist.append((orig.shape[i] - m1) / (dimensions[i] - m1) *
                       (base + offset) - offset)

    # specify old coordinates
    old_coords = [np.arange(i, dtype=np.float) for i in orig.shape]

    # first interpolation - for ndims = any
    mint = scipy.interpolate.interp1d(old_coords[-1], orig, bounds_error=False,
                                      fill_value=min(old_coords[-1]), kind=method)

    new_data = mint(dimlist[-1])

    trorder = [orig.ndim - 1] + list(range(orig.ndim - 1))
    for i in range(orig.ndim - 2, -1, -1):
        new_data = new_data.transpose(trorder)

        mint = scipy.interpolate.interp1d(old_coords[i], new_data,
            bounds_error=False, fill_value=min(old_coords[i]), kind=method)
        new_data = mint(dimlist[i])

    if orig.ndim > 1:
        # need one more transpose to return to original dimensions
        new_data = new_data.transpose(trorder)

    return new_data


def _resample_neighbor(orig, dimensions, offset, m1):
    """Resample Map using closest-value interpolation."""

    dimlist = []

    for i in range(orig.ndim):
        base = np.indices(dimensions)[i]
        dimlist.append((orig.shape[i] - m1) / (dimensions[i] - m1) *
                       (base + offset) - offset)
    cd = np.array(dimlist).round().astype(int)

    return orig[list(cd)]


def _resample_spline(orig, dimensions, offset, m1):
    """Resample Map using spline-based interpolation."""

    oslices = [slice(0, j) for j in orig.shape]
    # FIXME: not used?!
    old_coords = np.ogrid[oslices] #pylint: disable=W0612
    nslices = [slice(0, j) for j in list(dimensions)]
    newcoords = np.mgrid[nslices]

    newcoords_dims = list(range(np.rank(newcoords)))

    #make first index last
    newcoords_dims.append(newcoords_dims.pop(0))
    newcoords_tr = newcoords.transpose(newcoords_dims) #pylint: disable=W0612

    # makes a view that affects newcoords
    newcoords_tr += offset

    deltas = (np.asarray(orig.shape) - m1) / (dimensions - m1)
    newcoords_tr *= deltas

    newcoords_tr -= offset

    return scipy.ndimage.map_coordinates(orig, newcoords)


def reshape_image_to_4d_superpixel(img, dimensions, offset):
    """Re-shape the two dimension input image into a a four dimensional
    array whose first and third dimensions express the number of original
    pixels in the x and y directions that form one superpixel. The reshaping
    makes it very easy to perform operations on superpixels.

    An application of this reshaping is the following.  Let's say you have an
    array

    x = np.array([[0, 0, 0, 1, 1, 1],
                  [0, 0, 1, 1, 0, 0],
                  [1, 1, 0, 0, 1, 1],
                  [0, 0, 0, 0, 1, 1],
                  [1, 0, 1, 0, 1, 1],
                  [0, 0, 1, 0, 0, 0]])

    and you want to sum over 2x2 non-overlapping sub-arrays.  For example, you
    could have a noisy image and you want to increase the signal-to-noise ratio.
    Summing over all the non-overlapping 2x2 sub-arrays will create a
    superpixel array of the original data.  Every pixel in the superpixel array
    is the sum of the values in a 2x2 sub-array of the original array,

    This summing can be done by reshaping the array

    y = x.reshape(3,2,3,2)

    and then summing over the 1st and third directions

    y2 = y.sum(axis=3).sum(axis=1)

    which gives the expected array.

    array([[0, 3, 2],
           [2, 0, 4],
           [1, 2, 2]])

    Parameters
    ----------
    img : `numpy.ndarray`
        A two-dimensional `~numpy.ndarray` of the form (y, x).

    dimensions : array-like
        A two element array-like object containing integers that describe the
        superpixel summation in the (y, x) directions.

    offset : array-like
        A two element array-like object containing integers that describe
        where in the input image the array reshaping begins in the (y, x)
        directions.

    Returns
    -------
    A four dimensional `~numpy.ndarray` that can be used to easily create
    two-dimensional arrays of superpixels of the input image.

    References
    ----------
    Taken from
    http://mail.scipy.org/pipermail/numpy-discussion/2010-July/051760.html

    """
    # make sure the input dimensions are integers
    dimensions = [int(dim) for dim in dimensions]

    # New dimensions of the final image
    na = int(np.floor((img.shape[0] - offset[0]) / dimensions[0]))
    nb = int(np.floor((img.shape[1] - offset[1]) / dimensions[1]))

    # Reshape up to a higher dimensional array which is useful for higher
    # level operations
    return (img[offset[0]:offset[0] + na*dimensions[0],
                offset[1]:offset[1] + nb*dimensions[1]]).reshape(na, dimensions[0], nb, dimensions[1])


class UnrecognizedInterpolationMethod(ValueError):
    """Unrecognized interpolation method specified."""
    pass


class UnequalNumDimensions(ValueError):
    """Number of dimensions does not match input array"""
    pass
