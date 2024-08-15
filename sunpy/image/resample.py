"""
Image resampling methods.
"""
import numpy as np
import scipy.interpolate
import scipy.ndimage

__all__ = ['resample', 'reshape_image_to_4d_superpixel']


def resample(orig, dimensions, method='linear', center=False, minusone=False):
    """
    Returns a new `numpy.ndarray` that has been resampled up or down.

    Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same coordinate lookup points
    as IDL's ``congrid`` routine (which apparently originally came from a
    VAX/VMS routine of the same name.)

    Parameters
    ----------
    orig : `numpy.ndarray`
        Original input array.
    dimensions : `tuple`
        Dimensions that new `numpy.ndarray` should have.
    method : {``"nearest"``, ``"linear"``, ``"spline"``}, optional
        Method to use for resampling interpolation.
            * nearest and linear - Uses "n x 1D" interpolations calculated by
              `scipy.interpolate.interp1d`.
            * spline - Uses `scipy.ndimage.map_coordinates`
    center : `bool`, optional
        If `False` (default) the interpolation points are at the front edge of the bin.
        If `True`, interpolation points are at the centers of the bins
    minusone : `bool`, optional
        For ``orig.shape = (i,j)`` & new dimensions ``= (x,y)``, if set to `False`
        (default) ``orig`` is resampled by factors of ``(i/x) * (j/y)``,
        otherwise ``orig`` is resampled by ``(i-1)/(x-1) * (j-1)/(y-1)``.
        This prevents extrapolation one element beyond bounds of input array.

    Returns
    -------
    out : `numpy.ndarray`
        A new `numpy.ndarray` which has been resampled to the desired dimensions.

    References
    ----------
    https://scipy-cookbook.readthedocs.io/items/Rebinning.html
    """

    # Verify that number dimensions requested matches original shape
    if len(dimensions) != orig.ndim:
        raise UnequalNumDimensions("Number of dimensions must remain the same "
                                   "when calling resample.")

    # TODO: Will this be okay for integer (e.g. JPEG 2000) data?
    if orig.dtype not in [np.float64, np.float32]:
        orig = orig.astype(np.float64)

    dimensions = np.asarray(dimensions, dtype=np.float64)
    m1 = np.array(minusone, dtype=np.int64)  # array(0) or array(1)
    offset = np.float64(center * 0.5)       # float64(0.) or float64(0.5)

    # Resample data
    if method in ['nearest', 'linear']:
        data = _resample_nearest_linear(orig, dimensions, method,
                                        offset, m1)
    elif method == 'spline':
        data = _resample_spline(orig, dimensions, offset, m1)
    else:
        raise UnrecognizedInterpolationMethod(f"Unrecognized interpolation method requested: {method}")

    return data


def _resample_nearest_linear(orig, dimensions, method, offset, m1):
    """
    Resample Map using either linear or nearest interpolation.

    Parameters
    ----------
    orig : array-like
        Original data.
    dimensions : `tuple`
        Dimensions of resampled data.
    method : `str`
        Interpolation method passed to `~scipy.interpolate.interpn`
    offset : `float`
        Either 0 or 0.5, depending on whether interpolation is at the edge or
        centers of bins.
    m1 : 0 or 1
        For ``orig.shape = (i,j)`` & new dimensions ``= (x,y)``, if set to `False`
        (default) ``orig`` is resampled by factors of ``(i/x) * (j/y)``,
        otherwise ``orig`` is resampled by ``(i-1)/(x-1) * (j-1)/(y-1)``.
        This prevents extrapolation one element beyond bounds of input array.
    """
    old_coords = [np.arange(i, dtype=float) + offset for i in orig.shape]
    scale = (orig.shape - m1) / (dimensions - m1)
    new_coords = [(np.arange(dimensions[i], dtype=float) + offset) * scale[i] for i in
                  range(len(dimensions))]
    new_coords = np.stack(np.meshgrid(*new_coords, indexing='ij'), axis=-1)
    # fill_value = None extrapolates outside the domain
    new_data = scipy.interpolate.interpn(old_coords, orig, new_coords,
                                         method=method, bounds_error=False,
                                         fill_value=None)

    return new_data


def _resample_spline(orig, dimensions, offset, m1):
    """
    Resample Map using spline-based interpolation.
    """
    nslices = [slice(0, j) for j in list(dimensions)]
    newcoords = np.mgrid[nslices]

    newcoords_dims = list(range(newcoords.ndim))

    # make first index last
    newcoords_dims.append(newcoords_dims.pop(0))
    newcoords_tr = newcoords.transpose(newcoords_dims)

    # makes a view that affects newcoords
    newcoords_tr += offset

    deltas = (np.asarray(orig.shape) - m1) / (dimensions - m1)
    newcoords_tr *= deltas

    newcoords_tr -= offset

    return scipy.ndimage.map_coordinates(orig, newcoords)


def reshape_image_to_4d_superpixel(img, dimensions, offset):
    """
    Re-shape the two-dimensional input image into a four-dimensional array
    whose first and third dimensions express the number of original pixels in
    the "x" and "y" directions that form one superpixel. The reshaping makes it
    very easy to perform operations on superpixels.

    An application of this reshaping is the following. Let's say you have an
    array::

         x = np.array([[0, 0, 0, 1, 1, 1],
                       [0, 0, 1, 1, 0, 0],
                       [1, 1, 0, 0, 1, 1],
                       [0, 0, 0, 0, 1, 1],
                       [1, 0, 1, 0, 1, 1],
                       [0, 0, 1, 0, 0, 0]])

    and you want to sum over 2x2 non-overlapping sub-arrays. For example, you
    could have a noisy image and you want to increase the signal-to-noise ratio.
    Summing over all the non-overlapping 2x2 sub-arrays will create a
    superpixel array of the original data. Every pixel in the superpixel array
    is the sum of the values in a 2x2 sub-array of the original array.

    This summing can be done by reshaping the array::

         y = x.reshape(3,2,3,2)

    and then summing over the 1st and third directions::

         y2 = y.sum(axis=3).sum(axis=1)

    which gives the expected array::

         array([[0, 3, 2],
               [2, 0, 4],
               [1, 2, 2]])

    Parameters
    ----------
    img : `numpy.ndarray`
        A two-dimensional `numpy.ndarray` of the form ``(y, x)``.
    dimensions : array-like
        A two element array-like object containing integers that describe the
        superpixel summation in the ``(y, x)`` directions.
    offset : array-like
        A two element array-like object containing integers that describe
        where in the input image the array reshaping begins in the ``(y, x)``
        directions.

    Returns
    -------
    A four dimensional `numpy.ndarray` that can be used to easily create
    two-dimensional arrays of superpixels of the input image.
    """
    # make sure the input dimensions are integers
    dimensions = [int(dim) for dim in dimensions]

    # New dimensions of the final image
    na = int(np.floor((img.shape[0] - offset[0]) / dimensions[0]))
    nb = int(np.floor((img.shape[1] - offset[1]) / dimensions[1]))

    # Reshape up to a higher dimensional array which is useful for higher
    # level operations
    return (img[int(offset[0]):int(offset[0] + na * dimensions[0]),
                int(offset[1]):int(offset[1] + nb * dimensions[1])]).reshape(na, dimensions[0], nb, dimensions[1])


class UnrecognizedInterpolationMethod(ValueError):
    """
    Unrecognized interpolation method specified.
    """


class UnequalNumDimensions(ValueError):
    """
    Number of dimensions does not match input array.
    """
