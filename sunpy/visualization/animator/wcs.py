from functools import partial

import numpy as np

import astropy.units as u
from astropy.wcs.wcsapi import BaseLowLevelWCS

from sunpy.extern import modest_image
from sunpy.visualization.animator.base import ArrayAnimator


__all__ = ['ArrayAnimatorWCS']


class ArrayAnimatorWCS(ArrayAnimator):
    """
    Animate an array with associated `~astropy.wcs.wcsapi.BaseLowLevelWCS` object.

    The following keyboard shortcuts are defined in the viewer:

    * 'left': previous step on active slider.
    * 'right': next step on active slider.
    * 'top': change the active slider up one.
    * 'bottom': change the active slider down one.
    * 'p': play/pause active slider.

    Parameters
    ----------
    data: `numpy.ndarray`
        The data to be visualized.
    wcs: `astropy.wcs.wcsapi.BaseLowLevelWCS`
        The world coordinate object associated with the array.
    slices: `tuple` or `list`
        A list specifying which axes of the array should be plotted on which
        axes. The list should be the same length as the number of pixel
        dimensions with ``'x'`` and (optionally) ``'y'`` in the elements
        corresponding to the axes to be plotted. If only ``'x'`` is present a
        line plot will be drawn. All other elements should be `0`.
    coord_params: `dict`, optional
        This dict allows you to override
        `~astropy.visualization.wcsaxes.WCSAxes` parameters for each world
        coordinate. The keys of this dictionary should be a value which can be
        looked up in ``WCSAxes.coords`` (i.e. ``em.wl`` or ``hpln``) and the
        values should be a dict which supports the following keys, and passes
        their values to the associated `~astropy.visualization.wcsaxes.WCSAxes`
        methods.

        * ``format_unit``: `~astropy.visualization.wcsaxes.CoordinateHelper.set_format_unit`
        * ``major_formatter``: `~astropy.visualization.wcsaxes.CoordinateHelper.set_major_formatter`
        * ``axislabel``: `~astropy.visualization.wcsaxes.CoordinateHelper.set_axislabel`
        * ``grid``: `~astropy.visualization.wcsaxes.CoordinateHelper.grid` (The value should be a dict of keyword arguments to ``grid()`` or `True`).
        * ``ticks``: `dict` the keyword arguments to the `~astropy.visualization.wcsaxes.CoordinateHelper.set_ticks` method.
    ylim: `tuple` or `str`, optional
       The yaxis limits to use when drawing a line plot, if 'fixed' then use
       the global data limits, if 'dynamic' then set the y limit for each frame
       individually (meaning the y limits change as you animate).
    ylabel: `string`, optional
       The yaxis label to use when drawing a line plot. Setting the label on
       the y-axis on an image plot should be done via ``coord_params``.

    """

    def __init__(self, data, wcs, slices, coord_params=None, ylim='dynamic', ylabel=None, **kwargs):
        if not isinstance(wcs, BaseLowLevelWCS):
            raise ValueError("A WCS object should be provided that implements the astropy WCS API.")
        if wcs.pixel_n_dim != data.ndim:
            raise ValueError("Dimensionality of the data and WCS object do not match.")
        if len(slices) != wcs.pixel_n_dim:
            raise ValueError("slices should be the same length as the number of pixel dimensions.")
        if "x" not in slices:
            raise ValueError("slices should contain at least 'x' to indicate the axis to plot on the x axis.")

        self.plot_dimensionality = 1

        image_axes = [slices[::-1].index("x")]
        if "y" in slices:
            image_axes.append(slices[::-1].index("y"))
            self.plot_dimensionality = 2

        if self.plot_dimensionality == 1:
            try:
                from astropy.visualization.wcsaxes.frame import RectangularFrame1D
            except ImportError as e:
                raise ImportError("Astropy 4.0 must be installed to do line plotting with WCSAxes.") from e

        self.naxis = data.ndim
        self.num_sliders = self.naxis - self.plot_dimensionality
        self.slices_wcsaxes = list(slices)
        self.wcs = wcs
        self.coord_params = coord_params
        self.ylim = ylim
        self.ylabel = ylabel

        extra_slider_labels = []
        if "slider_functions" in kwargs and "slider_labels" not in kwargs:
            extra_slider_labels = [a.__name__ for a in kwargs['slider_functions']]


        slider_labels = [w for i, w in enumerate(self._get_wcs_labels()) if not slices[i]][::-1] + extra_slider_labels
        slider_labels = kwargs.pop("slider_labels", slider_labels)

        super().__init__(data, image_axes=image_axes, axis_ranges=None,
                         slider_labels=slider_labels,
                         **kwargs)


    def _get_wcs_labels(self):
        """
        Read first the axes names property of the wcs and fall back to physical types.
        """
        # world_axis_names was only added to the APE 14 API in 4.0, so do this for backwards compatibility.
        world_axis_names = getattr(self.wcs, "world_axis_names", [''] * self.wcs.world_n_dim)
        # Return the name if it is set, or the physical type if it is not.
        return [l or t for l, t in zip(world_axis_names, self.wcs.world_axis_physical_types)]

    def _partial_pixel_to_world(self, pixel_dimension, pixel_coord):
        """
        Return the world coordinate along one axis, if it is only
        correlated to that axis.
        """
        wcs_dimension = self.wcs.pixel_n_dim - pixel_dimension - 1
        corr = self.wcs.axis_correlation_matrix[:, wcs_dimension]

        # If more than one world axis is linked to this dimension we can't
        # display the world coordinate because we have no way of picking,
        # so we just display pixel index.
        if len(np.nonzero(corr)[0]) != 1:
            return pixel_coord * u.pix

        # We know that the coordinate we care about is independent of the
        # other axes, so we can set the pixel coordinates to 0.
        coords = [0] * self.wcs.pixel_n_dim
        coords[wcs_dimension] = pixel_coord
        wc = self.wcs.pixel_to_world_values(*coords)[wcs_dimension]
        return u.Quantity(wc, unit=self.wcs.world_axis_units[wcs_dimension])

    def _sanitize_axis_ranges(self, *args):
        """
        This overrides the behaviour of ArrayAnimator to generate axis_ranges
        based on the WCS.
        """

        axis_ranges = [None] * self.wcs.pixel_n_dim
        for i in self.slider_axes:
            axis_ranges[i] = partial(self._partial_pixel_to_world, i)

        return axis_ranges, None

    def _apply_coord_params(self, axes):
        if self.coord_params is None:
            return

        for coord_name in self.coord_params:
            coord = axes.coords[coord_name]
            params = self.coord_params[coord_name]

            format_unit = params.get("format_unit", None)
            if format_unit:
                coord.set_format_unit(format_unit)

            major_formatter = params.get("major_formatter", None)
            if major_formatter:
                coord.set_major_formatter(major_formatter)

            axislabel = params.get("axislabel", None)
            if axislabel:
                coord.set_axislabel(axislabel)

            grid = params.get("grid", None)
            if grid is not None:
                if not isinstance(grid, dict):
                    grid = {}
                coord.grid(**grid)

            ticks = params.get("ticks", None)
            if ticks is not None:
                if not isinstance(ticks, dict):
                    raise TypeError("The 'ticks' value in the coord_params dictionary must be a dict.")
                coord.set_ticks(**ticks)

    def _get_main_axes(self):
        axes = self.fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=self.wcs,
                                 slices=self.slices_wcsaxes)

        self._apply_coord_params(axes)

        return axes

    def plot_start_image(self, ax):
        if self.plot_dimensionality == 1:
            artist = self.plot_start_image_1d(ax)

        elif self.plot_dimensionality == 2:
            artist = self.plot_start_image_2d(ax)

        return artist

    def update_plot(self, val, artist, slider):
        """
        Update the plot when a slider changes.

        This method both updates the state of the Animator and also re-draws
        the matplotlib artist.
        """
        ind = int(val)
        ax_ind = self.slider_axes[slider.slider_ind]
        self.frame_slice[ax_ind] = ind
        self.slices_wcsaxes[self.wcs.pixel_n_dim - ax_ind - 1] = ind

        if self.plot_dimensionality == 1:
            self.update_plot_1d(val, artist, slider)
        elif self.plot_dimensionality == 2:
            self.update_plot_2d(val, artist, slider)

        self._apply_coord_params(self.axes)
        return super().update_plot(val, artist, slider)

    def plot_start_image_1d(self, ax):
        """
        Set up a line plot.

        When plotting with WCSAxes, we always plot against pixel coordinate.
        """
        if self.ylim != 'dynamic':
            ylim = self.ylim
            if ylim == 'fixed':
                ylim = (self.data.min(), self.data.max())
            ax.set_ylim(ylim)

        if self.ylabel:
            ax.set_ylabel(self.ylabel)

        line, = ax.plot(self.data[self.frame_index], **self.imshow_kwargs)
        return line

    @property
    def data_transposed(self):
        """
        Return data for 2D plotting, transposed if needed.
        """
        if self.slices_wcsaxes.index('y') < self.slices_wcsaxes.index("x"):
            return self.data[self.frame_index].transpose()
        else:
            return self.data[self.frame_index]

    def update_plot_1d(self, val, line, slider):
        """
        Update the line plot.
        """
        if val != slider.cval:
            self.axes.reset_wcs(wcs=self.wcs, slices=self.slices_wcsaxes)
            line.set_ydata(self.data[self.frame_index])

            # If we are not setting ylim globally then we set it per frame.
            if self.ylim == 'dynamic':
                self.axes.set_ylim(self.data[self.frame_index].min(),
                                   self.data[self.frame_index].max())
            slider.cval = val

    def plot_start_image_2d(self, ax):
        """
        Setup an image plot.
        """
        imshow_args = {'interpolation': 'nearest',
                       'origin': 'lower'}
        imshow_args.update(self.imshow_kwargs)

        im = modest_image.imshow(ax, self.data_transposed, **imshow_args)

        if 'extent' in imshow_args:
            ax.set_xlim(imshow_args['extent'][:2])
            ax.set_ylim(imshow_args['extent'][2:])
        else:
            ny, nx = self.data_transposed.shape
            ax.set_xlim(-0.5, nx - 0.5)
            ax.set_ylim(-0.5, ny - 0.5)

        ax.dataLim.intervalx = ax.get_xlim()
        ax.dataLim.intervaly = ax.get_ylim()

        if self.if_colorbar:
            self._add_colorbar(im)

        return im

    def update_plot_2d(self, val, im, slider):
        """
        Update the image plot.
        """
        if val != slider.cval:
            self.axes.reset_wcs(wcs=self.wcs, slices=self.slices_wcsaxes)
            im.set_array(self.data_transposed)
            slider.cval = val
