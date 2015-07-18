"""Image resampling methods"""
from __future__ import absolute_import

import numpy as np
import scipy.interpolate
import scipy.ndimage

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

    trorder = [orig.ndim - 1] + range(orig.ndim - 1)
    for i in xrange(orig.ndim - 2, -1, -1):
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

    for i in xrange(orig.ndim):
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

    newcoords_dims = range(np.rank(newcoords))

    #make first index last
    newcoords_dims.append(newcoords_dims.pop(0))
    newcoords_tr = newcoords.transpose(newcoords_dims) #pylint: disable=W0612

    # makes a view that affects newcoords
    newcoords_tr += offset

    deltas = (np.asarray(orig.shape) - m1) / (dimensions - m1)
    newcoords_tr *= deltas

    newcoords_tr -= offset

    return scipy.ndimage.map_coordinates(orig, newcoords)


def reshape_image_to_4d_superpixel(img, dimensions):
    """Re-shape the two dimension input image into a a four dimensional
    array whose 1st and third dimensions express the number of original
    pixels in the x and y directions form one superpixel. The reshaping
    makes it very easy to perform operations on superpixels.


    Parameters
    ----------
    img : `numpy.ndarray`
        A two-dimensional `~numpy.ndarray` of the form (y, x)

    dimensions : array-like
        A two element array-like object containing integers that describe the
        superpixel summation in the (y, x) directions

    Returns
    -------
    A four dimensional `~numpy.ndarray` that can be used to easily create
    two-dimensional arrays of superpixels of the input image

    References
    ----------
    Taken from
    http://mail.scipy.org/pipermail/numpy-discussion/2010-July/051760.html

    """
    # check that the dimensions divide into the image size exactly

    if np.all((np.array(img.shape) % np.array(dimensions) != 0)):
        raise ValueError('New dimensions must divide original image size exactly.')

    # Reshape up to a higher dimensional array which is useful for higher
    # level operations
    return img.reshape(img.shape[0] / dimensions[0],
                       dimensions[0],
                       img.shape[1] / dimensions[1],
                       dimensions[1])


class UnrecognizedInterpolationMethod(ValueError):
    """Unrecognized interpolation method specified."""
    pass


class UnequalNumDimensions(ValueError):
    """Number of dimensions does not match input array"""
    pass
