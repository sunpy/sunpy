import numpy as np

from astropy import units as u

from sunpy.image.resample import reshape_image_to_4d_superpixel
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

    @deprecated(since="6.1", message="sunpy.map.GenericMap.superpixel() will be removed in 6.1, use sunpy.map.GenericMap.rebin() instead")
    @u.quantity_input
    def superpixel(self, dimensions: u.pixel, offset: u.pixel = (0, 0)*u.pixel, func=np.sum):
        """Returns a new map consisting of superpixels formed by applying
        'func' to the original map data.

        Parameters
        ----------
        dimensions : tuple
            One superpixel in the new map is equal to (dimension[0],
            dimension[1]) pixels of the original map.
            The first argument corresponds to the 'x' axis and the second
            argument corresponds to the 'y' axis. If non-integer values are provided,
            they are rounded using `int`.
        offset : tuple
            Offset from (0,0) in original map pixels used to calculate where
            the data used to make the resulting superpixel map starts.
            If non-integer value are provided, they are rounded using `int`.
        func
            Function applied to the original data.
            The function 'func' must take a numpy array as its first argument,
            and support the axis keyword with the meaning of a numpy axis
            keyword (see the description of `~numpy.sum` for an example.)
            The default value of 'func' is `~numpy.sum`; using this causes
            superpixel to sum over (dimension[0], dimension[1]) pixels of the
            original map.

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new Map which has superpixels of the required size.

        References
        ----------
        | `Summarizing blocks of an array using a moving window <https://mail.scipy.org/pipermail/numpy-discussion/2010-July/051760.html>`_
        """

        # Note: because the underlying ndarray is transposed in sense when
        #   compared to the Map, the ndarray is transposed, resampled, then
        #   transposed back.
        # Note: "center" defaults to True in this function because data
        #   coordinates in a Map are at pixel centers.

        if (offset.value[0] < 0) or (offset.value[1] < 0):
            raise ValueError("Offset is strictly non-negative.")

        # These are rounded by int() in reshape_image_to_4d_superpixel,
        # so round here too for use in constructing metadata later.
        dimensions = [int(dim) for dim in dimensions.to_value(u.pix)]
        offset = [int(off) for off in offset.to_value(u.pix)]

        # Make a copy of the original data, perform reshaping, and apply the
        # function.
        if self.mask is not None:
            data = np.ma.array(self.data.copy(), mask=self.mask)
        else:
            data = self.data.copy()

        reshaped = reshape_image_to_4d_superpixel(data,
                                                  [dimensions[1], dimensions[0]],
                                                  [offset[1], offset[0]])
        new_array = func(func(reshaped, axis=3), axis=1)

        # Update image scale and number of pixels

        # create copy of new meta data
        new_meta = self.meta.copy()

        # Update metadata
        for key in {'cdelt1', 'cd1_1', 'cd2_1'} & self.meta.keys():
            new_meta[key] *= dimensions[0]
        for key in {'cdelt2', 'cd1_2', 'cd2_2'} & self.meta.keys():
            new_meta[key] *= dimensions[1]
        if 'pc1_1' in self.meta:
            new_meta['pc1_2'] *= dimensions[1] / dimensions[0]
            new_meta['pc2_1'] *= dimensions[0] / dimensions[1]

        new_meta['crpix1'] = ((self.reference_pixel.x.to_value(u.pix) +
                               0.5 - offset[0]) / dimensions[0]) + 0.5
        new_meta['crpix2'] = ((self.reference_pixel.y.to_value(u.pix) +
                               0.5 - offset[1]) / dimensions[1]) + 0.5
        new_meta['naxis1'] = new_array.shape[1]
        new_meta['naxis2'] = new_array.shape[0]

        # Create new map instance
        if self.mask is not None:
            new_data = np.ma.getdata(new_array)
            new_mask = np.ma.getmask(new_array)
        else:
            new_data = new_array
            new_mask = None

        # Create new map with the modified data
        new_map = self._new_instance(new_data, new_meta, self.plot_settings, mask=new_mask)
        return new_map
