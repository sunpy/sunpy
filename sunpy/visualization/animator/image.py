import matplotlib as mpl

from astropy.wcs.wcsapi import BaseLowLevelWCS

from sunpy.util.decorators import deprecated
from sunpy.visualization.animator.base import ArrayAnimator

__all__ = ['ImageAnimator', 'ImageAnimatorWCS']


class ImageAnimator(ArrayAnimator):
    """
    Create a matplotlib backend independent data explorer for 2D images.

    The following keyboard shortcuts are defined in the viewer:

    * 'left': previous step on active slider.
    * 'right': next step on active slider.
    * 'top': change the active slider up one.
    * 'bottom': change the active slider down one.
    * 'p': play/pause active slider.

    This viewer can have user defined buttons added by specifying the labels
    and functions called when those buttons are clicked as keyword arguments.

    Parameters
    ----------
    data: `numpy.ndarray`
        The data to be visualized.
    image_axes: `list`, optional
        A list of the axes order that make up the image.
    axis_ranges: `list` of physical coordinates for the `numpy.ndarray`, optional
        Defaults to `None` and array indices will be used for all axes.
        The `list` should contain one element for each axis of the `numpy.ndarray`.
        For the image axes a ``[min, max]`` pair should be specified which will be
        passed to `matplotlib.pyplot.imshow` as an extent.
        For the slider axes a ``[min, max]`` pair can be specified or an array the
        same length as the axis which will provide all values for that slider.

    Notes
    -----
    Extra keywords are passed to `~sunpy.visualization.animator.ArrayAnimator`.
    """
    def __init__(self, data, image_axes=[-2, -1], axis_ranges=None, **kwargs):
        # Check that number of axes is 2.
        if len(image_axes) != 2:
            raise ValueError("There can only be two spatial axes")
        # Define number of slider axes.
        self.naxis = data.ndim
        self.num_sliders = self.naxis-2
        # Define marker to determine if plot axes values are supplied via array of
        # pixel values or min max pair. This will determine the type of image produced
        # and hence how to plot and update it.
        self._non_regular_plot_axis = False
        # Run init for parent class
        super().__init__(data, image_axes=image_axes, axis_ranges=axis_ranges, **kwargs)

    def plot_start_image(self, ax):
        """
        Sets up plot of initial image.
        """
        # Create extent arg
        extent = []
        # reverse because numpy is in y-x and extent is x-y
        if max([len(self.axis_ranges[i]) for i in self.image_axes[::-1]]) > 2:
            self._non_regular_plot_axis = True
        for i in self.image_axes[::-1]:
            if self._non_regular_plot_axis is False and len(self.axis_ranges[i]) > 2:
                self._non_regular_plot_axis = True
            extent.append(self.axis_ranges[i][0])
            extent.append(self.axis_ranges[i][-1])

        imshow_args = {'interpolation': 'nearest',
                       'origin': 'lower'}
        imshow_args.update(self.imshow_kwargs)

        # If value along an axis is set with an array, generate a NonUniformImage
        if self._non_regular_plot_axis:
            # If user has inverted the axes, transpose the data so the dimensions match.
            if self.image_axes[0] < self.image_axes[1]:
                data = self.data[self.frame_index].transpose()
            else:
                data = self.data[self.frame_index]
            # Initialize a NonUniformImage with the relevant data and axis values and
            # add the image to the axes.
            im = mpl.image.NonUniformImage(ax, **imshow_args)
            im.set_data(self.axis_ranges[self.image_axes[0]], self.axis_ranges[self.image_axes[1]], data)
            ax.add_image(im)
            # Define the xlim and ylim from the pixel edges.
            ax.set_xlim(self.extent[0], self.extent[1])
            ax.set_ylim(self.extent[2], self.extent[3])
        else:
            # Else produce a more basic plot with regular axes.
            imshow_args.update({'extent': extent})
            im = ax.imshow(self.data[self.frame_index], **imshow_args)
        if self.if_colorbar:
            self._add_colorbar(im)

        return im

    def update_plot(self, val, im, slider):
        """
        Updates plot based on slider/array dimension being iterated.
        """
        ind = int(val)
        ax_ind = self.slider_axes[slider.slider_ind]
        self.frame_slice[ax_ind] = ind
        if val != slider.cval:
            if self._non_regular_plot_axis:
                if self.image_axes[0] < self.image_axes[1]:
                    data = self.data[self.frame_index].transpose()
                else:
                    data = self.data[self.frame_index]
                im.set_data(self.axis_ranges[self.image_axes[0]],
                            self.axis_ranges[self.image_axes[1]], data)
            else:
                im.set_array(self.data[self.frame_index])
            slider.cval = val
        # Update slider label to reflect real world values in axis_ranges.
        super().update_plot(val, im, slider)


@deprecated("1.1", alternative="sunpy.visualization.animator.ArrayAnimatorWCS")
class ImageAnimatorWCS(ImageAnimator):
    """
    Animates N-dimensional data with an associated World Coordinate System.

    The following keyboard shortcuts are defined in the viewer:

    * 'left': previous step on active slider.
    * 'right': next step on active slider.
    * 'top': change the active slider up one.
    * 'bottom': change the active slider down one.
    * 'p': play/pause active slider.

    This viewer can have user defined buttons added by specifying the labels
    and functions called when those buttons are clicked as keyword arguments.

    Parameters
    ----------
    data: `numpy.ndarray`
        The data to be visualized.
    wcs : `~astropy.wcs.wcsapi.BaseLowLevelWCS`
        The WCS object describing the physical coordinates of the data.
    image_axes: `list`, optional
        A list of the axes order that make up the image.
    unit_x_axis: `astropy.units.Unit`, optional
        The unit of X axis.
    unit_y_axis: `astropy.units.Unit`, optional
        The unit of Y axis.
    axis_ranges: `list`, optional
        Defaults to `None` and array indices will be used for all axes.
        The `list` should contain one element for each axis of the input data array.
        For the image axes a ``[min, max]`` pair should be specified which will be
        passed to `matplotlib.pyplot.imshow` as an extent.
        For the slider axes a ``[min, max]`` pair can be specified or an array the
        same length as the axis which will provide all values for that slider.

    Notes
    -----
    Extra keywords are passed to `~sunpy.visualization.animator.ArrayAnimator`.
    """
    def __init__(self, data, wcs, image_axes=[-1, -2], unit_x_axis=None, unit_y_axis=None,
                 axis_ranges=None, **kwargs):
        if not isinstance(wcs, BaseLowLevelWCS):
            raise ValueError("A WCS object should be provided that implements the astropy WCS API.")
        if wcs.pixel_n_dim is not data.ndim:
            raise ValueError("Dimensionality of the data and WCS object do not match.")
        self.wcs = wcs
        list_slices_wcsaxes = [0 for i in range(self.wcs.pixel_n_dim)]
        list_slices_wcsaxes[image_axes[0]] = 'x'
        list_slices_wcsaxes[image_axes[1]] = 'y'
        self.slices_wcsaxes = list_slices_wcsaxes[::-1]
        self.unit_x_axis = unit_x_axis
        self.unit_y_axis = unit_y_axis

        # Using `super()` here causes an error with the @deprecated decorator.
        ImageAnimator.__init__(self, data, image_axes=image_axes, axis_ranges=axis_ranges, **kwargs)

    def _get_main_axes(self):
        axes = self.fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=self.wcs,
                                 slices=self.slices_wcsaxes)
        self._set_unit_in_axis(axes)
        return axes

    def _set_unit_in_axis(self, axes):
        x_index = self.slices_wcsaxes.index("x")
        y_index = self.slices_wcsaxes.index("y")
        if self.unit_x_axis is not None:
            axes.coords[x_index].set_format_unit(self.unit_x_axis)
            axes.coords[x_index].set_ticks(exclude_overlapping=True)
        if self.unit_y_axis is not None:
            axes.coords[y_index].set_format_unit(self.unit_y_axis)
            axes.coords[y_index].set_ticks(exclude_overlapping=True)

    def plot_start_image(self, ax):
        """
        Sets up a plot of initial image.
        """
        imshow_args = {'interpolation': 'nearest',
                       'origin': 'lower',
                       }
        imshow_args.update(self.imshow_kwargs)
        im = ax.imshow(self.data[self.frame_index], **imshow_args)
        if self.if_colorbar:
            self._add_colorbar(im)
        return im

    def update_plot(self, val, im, slider):
        """
        Updates plot based on slider/array dimension being iterated.
        """
        ind = int(val)
        ax_ind = self.slider_axes[slider.slider_ind]
        self.frame_slice[ax_ind] = ind
        list_slices_wcsaxes = list(self.slices_wcsaxes)
        list_slices_wcsaxes[self.wcs.pixel_n_dim-ax_ind-1] = val
        self.slices_wcsaxes = list_slices_wcsaxes
        if val != slider.cval:
            self.axes.reset_wcs(wcs=self.wcs, slices=self.slices_wcsaxes)
            self._set_unit_in_axis(self.axes)
            im.set_array(self.data[self.frame_index])
            slider.cval = val
        # Update slider label to reflect real world values in axis_ranges.
        super().update_plot(val, im, slider)
