"""Image resampling methods"""
from __future__ import absolute_import

def resample(orig, dimensions, method='linear', center=False, minusone=False):
    """Returns a new ndarray that has been resampled up or down
    
    Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a 
    VAX/VMS routine of the same name.
    
    Parameters
    ----------
    dimensions : tuple
        Dimensions that new ndarray should have.
    method : {'neighbor' | 'nearest' | 'linear' | 'spline'}
        Method to use for resampling interpolation.
            * neighbor - Closest value from original data
            * nearest and linear - Uses n x 1-D interpolations using
              scipy.interpolate.interp1d
            * spline - Uses ndimage.map_coordinates
    center : bool
        If True, interpolation points are at the centers of the bins,
        otherwise points are at the front edge of the bin.
    minusone : bool
        For inarray.shape = (i,j) & new dimensions = (x,y), if set to False
        inarray is resampled by factors of (i/x) * (j/y), otherwise inarray 
        is resampled by(i-1)/(x-1) * (j-1)/(y-1)
        This prevents extrapolation one element beyond bounds of input 
        array.

    Returns
    -------
    out : ndarray
        A new ndarray which has been resampled to the desired dimensions.
    
    References
    ----------
    | http://www.scipy.org/Cookbook/Rebinning (Original source, 2011/11/19)
    """
    import numpy as np

    # Verify that number dimensions requested matches original shape
    if len(dimensions) != orig.ndim:
        raise UnequalNumDimensions("Number of dimensions must remain the same "
                                   "when calling resample.")

    #@note: will this be okay for integer (e.g. JPEG 2000) data?
    if not orig.dtype in [np.float64, np.float32]:
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
    """Resample Map using either linear or nearest interpolation"""
    import scipy.interpolate
    import numpy as np

    dimlist = []
    
    # calculate new dims
    for i in range(orig.ndim):
        base = np.arange(dimensions[i])
        dimlist.append((orig.shape[i] - m1) / (dimensions[i] - m1) *
                       (base + offset) - offset)

    # specify old coordinates
    old_coords = [np.arange(i, dtype=np.float) for i in orig.shape]

    # first interpolation - for ndims = any
    mint = scipy.interpolate.interp1d(old_coords[-1], orig, kind=method)
    new_data = mint(dimlist[-1])

    trorder = [orig.ndim - 1] + range(orig.ndim - 1)
    for i in xrange(orig.ndim - 2, -1, -1):
        new_data = new_data.transpose(trorder)

        mint = scipy.interpolate.interp1d(old_coords[i], new_data, 
                                          kind=method)
        new_data = mint(dimlist[i])

    if orig.ndim > 1:
        # need one more transpose to return to original dimensions
        new_data = new_data.transpose(trorder)

    return new_data

def _resample_neighbor(orig, dimensions, offset, m1):
    """Resample Map using closest-value interpolation"""
    import numpy as np
    
    dimlist = []
    
    for i in xrange(orig.ndim):
        base = np.indices(dimensions)[i]
        dimlist.append((orig.shape[i] - m1) / (dimensions[i] - m1) *
                       (base + offset) - offset)
    cd = np.array(dimlist).round().astype(int)
    
    return orig[list(cd)]

def _resample_spline(orig, dimensions, offset, m1):
    """Resample Map using spline-based interpolation"""
    import scipy.ndimage
    import numpy as np
    
    oslices = [slice(0, j) for j in orig.shape]
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


class UnrecognizedInterpolationMethod(ValueError):
    """Unrecognized interpolation method specified."""
    pass

class UnequalNumDimensions(ValueError):
    """Number of dimensions does not match input array"""
    pass
