import numpy as np

from . base import ArrayAnimator, edges_to_centers_nd

__all__ = ['LineAnimator']


class LineAnimator(ArrayAnimator):
    """
    Create a matplotlib backend independent data explorer for 1D plots.

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
        The y-axis data to be visualized

    plot_axis_index: `int`
        The axis used to plot against xdata.
        Default = -1, i.e. last dimension of arrary.

    fig: `matplotlib.figure`
        Figure to use

    axis_ranges: `list` of physical coordinates or None
        Each element of axis_ranges provides an array of physical coordinates for the
        pixel/bin edges along the corresponding dimension of the data array.
        If an element is None, array indices will be used for that axis.
        If axis_range itself is None, array indices will be used for all axes.
        For more information, see the Notes section of this docstring.

    xlabel: `str`
        Label of x-axis of plot.

    ylabel: `str`
        Label of y-axis of plot.

    xlim: `tuple`
        Limits of x-axis of plot.

    ylim: `tuple`
        Limits of y-axis of plot.

    interval: `int`
        Animation interval in ms

    button_labels: `list`
        List of strings to label buttons

    button_func: `list`
        List of functions to map to the buttons

    Extra keywords are passed to plot.

    Notes
    -----
    Additional information on API of axes_ranges kwarg.

    #. X-axis values must be supplied (if desired) as an array in the element of
       the axis_ranges list corresponding to the plot_axis_index in the data array,
       i.e. ``x_axis_values == axis_ranges[plot_axis_index]``

    #. The x-axis values represent the edges of the pixels/bins along the plotted
       axis, not the centers. Therefore there must be 1 more x-axis value than
       there are data points along the x-axis.

    #. The shape of the x-axis values array can take two forms.

       a) First, it can have a length 1 greater than the length of the data array
          along the dimension corresponding to the x-axis, i.e.
          ``len(axis_ranges[plot_axis_index]) == len(data[plot_axis_index])+1``
          In this scenario the same x-axis values are used in every frame of the
          animation.
       b) Second, the x-axis array can have the same shape as the data array, with
          the exception of the plotted axis which, as above, must be 1 greater than
          the length of the data array along that dimension.
          In this scenario the x-axis is refreshed for each frame. For example, if
          ``data.shape == axis_ranges[plot_axis_index] == (4, 3)``,
          where ``plot_axis_index == 0``, the 0th frame of the animation will show data from
          ``data[:, 0]`` with the x-axis described by ``axis_ranges[plot_axis_index][:, 0]``,
          while the 1st frame will show data from ``data[:, 1]`` with the x-axis described by
          ``axis_ranges[plot_axis_index][:, 1]``.

    #. The API holds for slider axes.

    """

    def __init__(self, data, plot_axis_index=-1, axis_ranges=None, ylabel=None, xlabel=None,
                 xlim=None, ylim=None, **kwargs):
        # Check inputs.
        self.plot_axis_index = int(plot_axis_index)
        if self.plot_axis_index not in range(-data.ndim, data.ndim):
            raise ValueError("plot_axis_index must be within range of number of data dimensions"
                             " (or equivalent negative indices).")
        if data.ndim < 2:
            raise ValueError("data must have at least two dimensions.  One for data "
                             "for each single plot and at least one for time/iteration.")
        # Define number of slider axes.
        self.naxis = data.ndim
        self.num_sliders = self.naxis-1
        # Attach data to class.
        if axis_ranges is not None and all(axis_range is None for axis_range in axis_ranges):
            axis_ranges = None
        if axis_ranges is None or axis_ranges[self.plot_axis_index] is None:
            self.xdata = np.arange(data.shape[self.plot_axis_index])
        else:
            # Else derive the xdata as the centers of the pixel/bin edges
            # supplied by the user for the plotted axis.
            self.xdata = edges_to_centers_nd(np.asarray(axis_ranges[self.plot_axis_index]),
                                             plot_axis_index)
        if ylim is None:
            ylim = (data.min(), data.max())
        if xlim is None:
            xlim = (self.xdata.min(), self.xdata.max())
        self.ylim = ylim
        self.xlim = xlim
        self.xlabel = xlabel
        self.ylabel = ylabel
        # Run init for base class
        super().__init__(data, image_axes=[self.plot_axis_index], axis_ranges=axis_ranges,
                         **kwargs)

    def plot_start_image(self, ax):
        """Sets up plot of initial image."""
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        if self.xlabel is not None:
            ax.set_xlabel(self.xlabel)
        if self.ylabel is not None:
            ax.set_ylabel(self.ylabel)
        plot_args = {}
        plot_args.update(self.imshow_kwargs)
        if self.xdata.shape == self.data.shape:
            item = [0] * self.data.ndim
            item[self.plot_axis_index] = slice(None)
            xdata = np.squeeze(self.xdata[tuple(item)])
        else:
            xdata = self.xdata
        line, = ax.plot(xdata, self.data[self.frame_index], **plot_args)
        return line

    def update_plot(self, val, line, slider):
        """Updates plot based on slider/array dimension being iterated."""
        val = int(val)
        ax_ind = self.slider_axes[slider.slider_ind]
        ind = int(np.argmin(np.abs(self.axis_ranges[ax_ind] - val)))
        self.frame_slice[ax_ind] = ind
        if val != slider.cval:
            line.set_ydata(self.data[self.frame_index])
            if self.xdata.shape == self.data.shape:
                item = [int(slid._slider.val) for slid in self.sliders]
                item[ax_ind] = val
                if self.plot_axis_index < 0:
                    i = self.data.ndim + self.plot_axis_index
                else:
                    i = self.plot_axis_index
                item.insert(i, slice(None))
                line.set_xdata(self.xdata[tuple(item)])
            slider.cval = val
        # Update slider label to reflect real world values in axis_ranges.
        super().update_plot(val, slider)
