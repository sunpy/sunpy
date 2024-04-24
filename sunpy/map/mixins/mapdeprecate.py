import numpy as np
from mapmeta import PixelPair

from astropy import units as u

__all__ = ['MapDeprecateMixin']

class MapDeprecateMixin:
    """
    This class contains the deprecated attributes on map
    """

    # Some numpy extraction
    @property
    def dimensions(self):
        """
        The dimensions of the array (x axis first, y axis second).
        """
        return PixelPair(*u.Quantity(np.flipud(self.data.shape), 'pixel'))

    @property
    def dtype(self):
        """
        The `numpy.dtype` of the array of the map.
        """
        return self.data.dtype

    @property
    def ndim(self):
        """
        The value of `numpy.ndarray.ndim` of the data array of the map.
        """
        return self.data.ndim

    def std(self, *args, **kwargs):
        """
        Calculate the standard deviation of the data array, ignoring NaNs.
        """
        return np.nanstd(self.data, *args, **kwargs)

    def mean(self, *args, **kwargs):
        """
        Calculate the mean of the data array, ignoring NaNs.
        """
        return np.nanmean(self.data, *args, **kwargs)

    def min(self, *args, **kwargs):
        """
        Calculate the minimum value of the data array, ignoring NaNs.
        """
        return np.nanmin(self.data, *args, **kwargs)

    def max(self, *args, **kwargs):
        """
        Calculate the maximum value of the data array, ignoring NaNs.
        """
        return np.nanmax(self.data, *args, **kwargs)
