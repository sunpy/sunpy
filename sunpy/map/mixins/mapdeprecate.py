import numpy as np

from astropy import units as u

from sunpy.util.decorators import deprecated
from .mapmeta import PixelPair

__all__ = ['MapDeprecateMixin']

class MapDeprecateMixin:
    """
    This class contains the deprecated attributes on map.
    """

    @property
    @deprecated(since="6.1", message=".dimensions will be removed in 6.1. Use sunpy.map.GenericMap.shape instead")
    def dimensions(self):
        """
        The dimensions of the array (x axis first, y axis second).
        """
        return PixelPair(*u.Quantity(np.flipud(self.data.shape), 'pixel'))


    @property
    @deprecated(since="6.1", message="sunpy.map.GenericMap.dtype will be removed in 6.1, use sunpy.map.GenericMap.data.dtype instead")
    def dtype(self):
        """
        The `numpy.dtype` of the array of the map.
        """
        return self.data.dtype

    @property
    @deprecated(since="6.1", message="sunpy.map.GenericMap.ndim will be removed in 6.1, use sunpy.map.GenericMap.data.ndim() instead")
    def ndim(self):
        """
        The value of `numpy.ndarray.ndim` of the data array of the map.
        """
        return self.data.ndim


    @deprecated(since="6.1", message="sunpy.map.GenericMap.std() will be removed in 6.1, use sunpy.map.GenericMap.data.std() instead")
    def std(self, *args, **kwargs):
        """
        Calculate the standard deviation of the data array, ignoring NaNs.
        """
        return np.nanstd(self.data, *args, **kwargs)

    @deprecated(since="6.1", message="sunpy.map.GenericMap.mean() will be removed in 6.1, use sunpy.map.GenericMap.data.mean() instead")
    def mean(self, *args, **kwargs):
        """
        Calculate the mean of the data array, ignoring NaNs.
        """
        return np.nanmean(self.data, *args, **kwargs)

    @deprecated(since="6.1", message="sunpy.map.GenericMap.min() will be removed in 6.1, use sunpy.map.GenericMap.data.min() instead")
    def min(self, *args, **kwargs):
        """
        Calculate the minimum value of the data array, ignoring NaNs.
        """
        return np.nanmin(self.data, *args, **kwargs)

    @deprecated(since="6.1", message="sunpy.map.GenericMap.max() will be removed in 6.1, use sunpy.map.GenericMap.data.max() instead")
    def max(self, *args, **kwargs):
        """
        Calculate the maximum value of the data array, ignoring NaNs.
        """
        return np.nanmax(self.data, *args, **kwargs)
