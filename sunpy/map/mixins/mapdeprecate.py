import numpy as np

from astropy import units as u

from sunpy.util.decorators import deprecated
from .mapmeta import PixelPair

__all__ = ['MapDeprecateMixin']

class MapDeprecateMixin:
    """
    This class contains the deprecated attributes on map
    """

    # Some numpy extraction
    @deprecated(since="6.1", message="this will be removed.", alternative="sunpy.map.GenericMap.shape")
    @property
    def dimensions(self):
        """
        The dimensions of the array (x axis first, y axis second).
        """
        return PixelPair(*u.Quantity(np.flipud(self.data.shape), 'pixel'))


    @deprecated(since="6.1", message="this will be removed.", alternative="sunpy.map.GenericMap.data.dtype")
    @property
    def dtype(self):
        """
        The `numpy.dtype` of the array of the map.
        """
        return self.data.dtype

    @deprecated(since="6.1", message="this will be removed.", alternative="sunpy.map.GenericMap.data.ndim()")
    @property
    def ndim(self):
        """
        The value of `numpy.ndarray.ndim` of the data array of the map.
        """
        return self.data.ndim

    @deprecated(since="6.1", message="this will be removed.", alternative="sunpy.map.GenericMap.data.std()")
    def std(self, *args, **kwargs):
        """
        Calculate the standard deviation of the data array, ignoring NaNs.
        """
        return np.nanstd(self.data, *args, **kwargs)

    @deprecated(since="6.1", message="this will be removed.", alternative="sunpy.map.GenericMap.data.mean()")
    def mean(self, *args, **kwargs):
        """
        Calculate the mean of the data array, ignoring NaNs.
        """
        return np.nanmean(self.data, *args, **kwargs)

    @deprecated(since="6.1", message="this will be removed.", alternative="sunpy.map.GenericMap.data.min()")
    def min(self, *args, **kwargs):
        """
        Calculate the minimum value of the data array, ignoring NaNs.
        """
        return np.nanmin(self.data, *args, **kwargs)

    @deprecated(since="6.1", message="this will be removed.", alternative="sunpy.map.GenericMap.data.max()")
    def max(self, *args, **kwargs):
        """
        Calculate the maximum value of the data array, ignoring NaNs.
        """
        return np.nanmax(self.data, *args, **kwargs)
