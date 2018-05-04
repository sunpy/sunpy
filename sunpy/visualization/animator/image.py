import numpy as np
import matplotlib as mpl
import astropy.wcs

from . base import ArrayAnimator

__all__ = ['ImageAnimator', 'ImageAnimatorWCS']


class ImageAnimator(ArrayAnimator):
    """
    Create a matplotlib backend independent data explorer for 2D images.

    The following keyboard shortcuts are defined in the viewer:

    - 'left': previous step on active slider
    - 'right': next step on active slider
    - 'top': change the active slider up one
    - 'bottom': change the active slider down one
    - 'p': play/pause active slider

    This viewer can have user defined buttons added by specifying the labels
    and functions called when those buttons are clicked as keyword arguments.

    Parameters
    ----------
    data: ndarray
        The data to be visualized >= 2D

    image_axes: list
        The two axes that make the image

    fig: mpl.figure
        Figure to use

    axis_ranges: list of physical coordinates for array or None
        If None array indices will be used for all axes.
        If a list it should contain one element for each axis of the numpy array.
        For the image axes a [min, max] pair should be specified which will be
        passed to :func:`matplotlib.pyplot.imshow` as extent.
        For the slider axes a [min, max] pair can be specified or an array the
        same length as the axis which will provide all values for that slider.
        If None is specified for an axis then the array indices will be used
        for that axis.

    interval: int
        Animation interval in ms

    colorbar: bool
        Plot colorbar

    button_labels: list
        List of strings to label buttons

    button_func: list
        List of functions to map to the buttons

    Extra keywords are passed to imshow.

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
        super(ImageAnimator, self).__init__(data, image_axes=image_axes,
                                            axis_ranges=axis_ranges, **kwargs)

    def plot_start_image(self, ax):
        """Sets up plot of initial image."""
        # Create extent arg
        extent = []
        # reverse because numpy is in y-x and extent is x-y
        for i in self.image_axes[::-1]:
            if self._non_regular_plot_axis is False and len(self.axis_ranges[i]) > 2:
                self._non_regular_plot_axis = True
            extent.append(self.axis_ranges[i][0])
            extent.append(self.axis_ranges[i][-1])

        imshow_args = {'interpolation': 'nearest',
                       'origin': 'lower',
                       'extent': extent}
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
            im.set_data(self.axis_ranges[self.image_axes[0]],
                        self.axis_ranges[self.image_axes[1]], data)
            ax.add_image(im)
            # Define the xlim and ylim bearing in mind that the pixel values along
            # the axes are plotted in the middle of the pixel.  Therefore make sure
            # there's half a pixel buffer either side of the ends of the axis ranges.
            x_ranges = self.axis_ranges[self.image_axes[0]]
            xlim = (x_ranges[0] - (x_ranges[1] - x_ranges[0]) / 2.,
                    x_ranges[-1] + (x_ranges[-1] - x_ranges[-2]) / 2.)
            y_ranges = self.axis_ranges[self.image_axes[1]]
            ylim = (y_ranges[0] - (y_ranges[1] - y_ranges[0]) / 2.,
                    y_ranges[-1] + (y_ranges[-1] - y_ranges[-2]) / 2.)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        else:
            # Else produce a more basic plot with regular axes.
            im = ax.imshow(self.data[self.frame_index], **imshow_args)
        if self.if_colorbar:
            self._add_colorbar(im)

        return im

    def update_plot(self, val, im, slider):
        """Updates plot based on slider/array dimension being iterated."""
        val = int(val)
        ax_ind = self.slider_axes[slider.slider_ind]
        ind = np.argmin(np.abs(self.axis_ranges[ax_ind] - val))
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


class ImageAnimatorWCS(ImageAnimator):
    """
    Animates N-dimensional data with the associated astropy WCS object.

    The following keyboard shortcuts are defined in the viewer:

    - 'left': previous step on active slider
    - 'right': next step on active slider
    - 'top': change the active slider up one
    - 'bottom': change the active slider down one
    - 'p': play/pause active slider

    This viewer can have user defined buttons added by specifying the labels
    and functions called when those buttons are clicked as keyword arguments.

    Parameters
    ----------
    data: `numpy.ndarray`
        The data to be visualized >= 2D

    wcs: `astropy.wcs.WCS`
        The wcs data.

    image_axes: `list`
        The two axes that make the image

    fig: `matplotlib.figure.Figure`
        Figure to use

    axis_ranges: list of physical coordinates for array or None
        If None array indices will be used for all axes.
        If a list it should contain one element for each axis of the numpy array.
        For the image axes a [min, max] pair should be specified which will be
        passed to :func:`matplotlib.pyplot.imshow` as extent.
        For the slider axes a [min, max] pair can be specified or an array the
        same length as the axis which will provide all values for that slider.
        If None is specified for an axis then the array indices will be used
        for that axis.

    interval: `int`
        Animation interval in ms

    colorbar: `bool`
        Plot colorbar

    button_labels: `list`
        List of strings to label buttons

    button_func: `list`
        List of functions to map to the buttons

    unit_x_axis: `astropy.units.Unit`
        The unit of x axis.

    unit_y_axis: `astropy.units.Unit`
        The unit of y axis.

    Extra keywords are passed to imshow.

    """
    def __init__(self, data, wcs=None, image_axes=[-1, -2], unit_x_axis=None, unit_y_axis=None,
                 axis_ranges=None, **kwargs):
        if not isinstance(wcs, astropy.wcs.WCS):
            raise ValueError("wcs data should be provided.")
        if wcs.wcs.naxis is not data.ndim:
            raise ValueError("Dimensions of data and wcs not matching")
        self.wcs = wcs
        list_slices_wcsaxes = [0 for i in range(self.wcs.naxis)]
        list_slices_wcsaxes[image_axes[0]] = 'x'
        list_slices_wcsaxes[image_axes[1]] = 'y'
        self.slices_wcsaxes = list_slices_wcsaxes[::-1]
        self.unit_x_axis = unit_x_axis
        self.unit_y_axis = unit_y_axis
        super(ImageAnimatorWCS, self).__init__(data, image_axes=image_axes,
                                               axis_ranges=axis_ranges, **kwargs)

    def _get_main_axes(self):
        axes = self.fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=self.wcs,
                                 slices=self.slices_wcsaxes)
        self._set_unit_in_axis(axes)
        return axes

    def _set_unit_in_axis(self, axes):
        if self.unit_x_axis is not None:
            axes.coords[2].set_format_unit(self.unit_x_axis)
            axes.coords[2].set_ticks(exclude_overlapping=True)
        if self.unit_y_axis is not None:
            axes.coords[1].set_format_unit(self.unit_y_axis)
            axes.coords[1].set_ticks(exclude_overlapping=True)

    def plot_start_image(self, ax):
        """Sets up plot of initial image."""
        imshow_args = {'interpolation': 'nearest',
                       'origin': 'lower',
                       }
        imshow_args.update(self.imshow_kwargs)
        im = ax.imshow(self.data[self.frame_index], **imshow_args)
        if self.if_colorbar:
            self._add_colorbar(im)
        return im

    def update_plot(self, val, im, slider):
        """Updates plot based on slider/array dimension being iterated."""
        val = int(val)
        ax_ind = self.slider_axes[slider.slider_ind]
        ind = np.argmin(np.abs(self.axis_ranges[ax_ind] - val))
        self.frame_slice[ax_ind] = ind
        list_slices_wcsaxes = list(self.slices_wcsaxes)
        list_slices_wcsaxes[self.wcs.naxis-ax_ind-1] = val
        self.slices_wcsaxes = list_slices_wcsaxes
        if val != slider.cval:
            self.axes.reset_wcs(wcs=self.wcs, slices=self.slices_wcsaxes)
            self._set_unit_in_axis(self.axes)
            im.set_array(self.data[self.frame_index])
            slider.cval = val
