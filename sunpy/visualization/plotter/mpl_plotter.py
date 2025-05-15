"""
Matplotlib plotter class
"""
import copy

import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import Longitude
from astropy.visualization import AsymmetricPercentileInterval
from astropy.visualization.wcsaxes import Quadrangle

import sunpy.visualization.drawing
from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.util.exceptions import warn_user
from sunpy.visualization import axis_labels_from_ctype, peek_show, wcsaxes_compat
from sunpy.visualization.colormaps import cm as sunpy_cm

__all__ = ['MapPlotter']

class MapPlotter:

    def __init__(self, smap, plot_settings=None):
        self.smap = smap
        self.plot_settings = plot_settings
        # Try and set the colormap. This is not always possible if this method
        # is run before map sources fix some of their metadata, so
        # just ignore any exceptions raised.
        try:
            cmap = self.smap._get_cmap_name()
            if cmap in sunpy_cm.cmlist:
                self.plot_settings['cmap'] = cmap
        except Exception:
            pass

    def _set_symmetric_vmin_vmax(self):
        """
        Set symmetric vmin and vmax about zero
        """
        threshold = np.nanmax(abs(self.smap.data))
        self.plot_settings['norm'].vmin = -threshold
        self.plot_settings['norm'].vmax = threshold

    @property
    def plot_settings(self):
        return self._plot_settings

    @plot_settings.setter
    def plot_settings(self, value):
        if self.smap.data.dtype == np.uint8:
            norm = None
        else:
            # Put import here to reduce sunpy.map import time
            from matplotlib import colors
            norm = colors.Normalize()
        # Visualization attributes
        self._plot_settings = {
            'cmap': 'gray',
            'norm': norm,
            'interpolation': 'nearest',
            'origin': 'lower',
        }
        if value:
            self._plot_settings.update(value)

    @property
    def cmap(self):
        """
        Return the `matplotlib.colors.Colormap` instance this map uses.
        """
        cmap = self.plot_settings['cmap']
        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)
            # Set the colormap to be this specific instance so we are not
            # returning a copy
            self.plot_settings['cmap'] = cmap
        return cmap

    @u.quantity_input
    def draw_grid(self, axes=None, grid_spacing: u.deg = 15*u.deg, annotate=True,
                  system='stonyhurst', **kwargs):
        """
        Draws a coordinate overlay on the plot in the Heliographic Stonyhurst
        coordinate system.

        To overlay other coordinate systems see the `WCSAxes Documentation
        <https://docs.astropy.org/en/stable/visualization/wcsaxes/overlaying_coordinate_systems.html>`_

        Parameters
        ----------
        axes : `~matplotlib.axes` or `None`
            Axes to plot limb on, or `None` to use current axes.
        grid_spacing : `~astropy.units.Quantity`
            Spacing for longitude and latitude grid, if length two it specifies
            (lon, lat) spacing.
        annotate : `bool`
            Passing `False` disables the axes labels and the ticks on the top and right axes.
        system : str
            Coordinate system for the grid. Must be 'stonyhurst' or 'carrington'.
        kwargs :
            Additional keyword arguments are passed to `~sunpy.visualization.wcsaxes_compat.wcsaxes_heliographic_overlay`.

        Returns
        -------
        overlay: `~astropy.visualization.wcsaxes.CoordinatesMap`
            The wcsaxes coordinate overlay instance.

        Notes
        -----
        Keyword arguments are passed onto the `sunpy.visualization.wcsaxes_compat.wcsaxes_heliographic_overlay` function.
        """
        axes = self._check_axes(axes)
        return wcsaxes_compat.wcsaxes_heliographic_overlay(axes,
                                                           grid_spacing=grid_spacing,
                                                           annotate=annotate,
                                                           obstime=self.smap.date,
                                                           rsun=self.smap.rsun_meters,
                                                           observer=self.smap.observer_coordinate,
                                                           system=system,
                                                           **kwargs)

    def draw_limb(self, axes=None, *, resolution=1000, **kwargs):
        """
        Draws the solar limb as seen by the map's observer.

        The limb is a circle for only the simplest plots.  If the coordinate frame of
        the limb is different from the coordinate frame of the plot axes, not only
        may the limb not be a true circle, a portion of the limb may be hidden from
        the observer.  In that case, the circle is divided into visible and hidden
        segments, represented by solid and dotted lines, respectively.

        Parameters
        ----------
        axes : `~matplotlib.axes` or ``None``
            Axes to plot limb on or ``None`` to use current axes.
        resolution : `int`
            The number of points to use to represent the limb.

        Returns
        -------
        visible : `~matplotlib.patches.Polygon` or `~matplotlib.patches.Circle`
            The patch added to the axes for the visible part of the limb (i.e., the
            "near" side of the Sun).
        hidden : `~matplotlib.patches.Polygon` or None
            The patch added to the axes for the hidden part of the limb (i.e., the
            "far" side of the Sun).

        Notes
        -----
        Keyword arguments are passed onto the patches.

        If the limb is a true circle, ``visible`` will instead be
        `~matplotlib.patches.Circle` and ``hidden`` will be ``None``. If there
        are no visible points (e.g., on a synoptic map any limb is fully)
        visible ``hidden`` will be ``None``.

        To avoid triggering Matplotlib auto-scaling, these patches are added as
        artists instead of patches.  One consequence is that the plot legend is not
        populated automatically when the limb is specified with a text label.  See
        :ref:`sphx_glr_gallery_text_labels_and_annotations_custom_legends.py` in
        the Matplotlib documentation for examples of creating a custom legend.
        """
        axes = self._check_axes(axes)
        return sunpy.visualization.drawing.limb(
            axes,
            self.smap.observer_coordinate,
            resolution=resolution,
            rsun=self.smap.rsun_meters,
            **kwargs
        )

    @u.quantity_input
    def draw_quadrangle(self, bottom_left, *, width: (u.deg, u.pix) = None, height: (u.deg, u.pix) = None,
                        axes=None, top_right=None, **kwargs):
        """
        Draw a quadrangle defined in world coordinates on the plot using Astropy's
        `~astropy.visualization.wcsaxes.Quadrangle`.

        This draws a quadrangle that has corners at ``(bottom_left, top_right)``,
        and has sides aligned with the coordinate axes of the frame of ``bottom_left``,
        which may be different from the coordinate axes of the map.

        If ``width`` and ``height`` are specified, they are respectively added to the
        longitude and latitude of the ``bottom_left`` coordinate to calculate a
        ``top_right`` coordinate.

        Parameters
        ----------
        bottom_left : `~astropy.coordinates.SkyCoord` or `~astropy.units.Quantity`
            The bottom-left coordinate of the rectangle. If a `~astropy.coordinates.SkyCoord` it can
            have shape ``(2,)`` and simultaneously define ``top_right``. If specifying
            pixel coordinates it must be given as an `~astropy.units.Quantity`
            object with pixel units (e.g., ``pix``).
        top_right : `~astropy.coordinates.SkyCoord` or `~astropy.units.Quantity`, optional
            The top-right coordinate of the quadrangle. If ``top_right`` is
            specified ``width`` and ``height`` must be omitted.
        width : `astropy.units.Quantity`, optional
            The width of the quadrangle. Required if ``top_right`` is omitted.
        height : `astropy.units.Quantity`
            The height of the quadrangle. Required if ``top_right`` is omitted.
        axes : `matplotlib.axes.Axes`
            The axes on which to plot the quadrangle. Defaults to the current
            axes.

        Returns
        -------
        quad : `~astropy.visualization.wcsaxes.Quadrangle`
            The added patch.

        Notes
        -----
        Extra keyword arguments to this function are passed through to the
        `~astropy.visualization.wcsaxes.Quadrangle` instance.

        Examples
        --------
        .. minigallery:: sunpy.map.GenericMap.draw_quadrangle
        """
        axes = self._check_axes(axes)

        if isinstance(bottom_left, u.Quantity):
            anchor, _, top_right, _ = self.smap._parse_submap_quantity_input(bottom_left, top_right, width, height)
            width, height = top_right - anchor
            transform = axes.get_transform(self.smap.wcs if self.smap.wcs is not axes.wcs else 'pixel')
            kwargs.update({"vertex_unit": u.pix})

        else:
            bottom_left, top_right = get_rectangle_coordinates(
                bottom_left, top_right=top_right, width=width, height=height)

            width = Longitude(top_right.spherical.lon - bottom_left.spherical.lon)
            height = top_right.spherical.lat - bottom_left.spherical.lat
            anchor = self.smap._get_lon_lat(bottom_left)
            transform = axes.get_transform(bottom_left.replicate_without_data())

        kwergs = {
            "transform": transform,
            "edgecolor": "white",
            "fill": False,
        }
        kwergs.update(kwargs)
        quad = Quadrangle(anchor, width, height, **kwergs)
        axes.add_patch(quad)
        return quad

    def _update_contour_args(self, contour_args):
        """
        Updates ``contour_args`` with values from ``plot_settings``.

        Parameters
        ----------
        contour_args : dict
            A dictionary of arguments to be used for contour plotting.

        Returns
        -------
        dict
            The updated ``contour_args`` dictionary.

        Notes
        -----
        - 'cmap': Set to `None` to avoid the error "ValueError: Either colors or cmap must be None".
        - 'interpolation': Removed because Matplotlib's contour function raises the warning
        "The following kwargs were not used by contour: 'interpolation'".
        - 'origin': If `'origin': 'lower'` is present, it is replaced with `'origin': None`,
        as `None` is the default value for Matplotlib's contour plots.
        """
        plot_settings = self.plot_settings.copy()
        contour_args_copy = contour_args.copy()
        contour_args.update(plot_settings)
        # Define default settings for normal plots and contour-specific updates
        original_plot_defaults = {
            'origin': 'lower',
        }
        default_contour_param = {
            'origin': None,
        }
        # Replace conflicting settings with contour defaults
        for key in original_plot_defaults:
            if key in contour_args and contour_args[key] == original_plot_defaults[key]:
                contour_args[key] = default_contour_param[key]
        # 'cmap' cannot be used for contour plots when levels are not None,
        # which is the case in composite maps.
        contour_args['cmap'] = None
        # custom 'norm' cannot be passed through plot_settings
        contour_args['norm'] = None
        # If 'draw_contour' is used, setting 'norm' and 'cmap' to None ensures the method arguments are applied.
        contour_args.update(contour_args_copy)
        contour_args.pop('interpolation')
        return contour_args


    def draw_contours(self, levels, axes=None, *, fill=False, **contour_args):
        """
        Draw contours of the data.

        Parameters
        ----------
        levels : `~astropy.units.Quantity`
            A list of numbers indicating the contours to draw. These are given
            as a percentage of the maximum value of the map data, or in units
            equivalent to the `~sunpy.map.GenericMap.unit` attribute.
        axes : `matplotlib.axes.Axes`
            The axes on which to plot the contours. Defaults to the current
            axes.
        fill : `bool`, optional
            Determines the style of the contours:
            - If `False` (default), contours are drawn as lines using :meth:`~matplotlib.axes.Axes.contour`.
            - If `True`, contours are drawn as filled regions using :meth:`~matplotlib.axes.Axes.contourf`.

        Returns
        -------
        cs : `list`
            The `~matplotlib.contour.QuadContourSet` object, after it has been added to
            ``axes``.

        Notes
        -----
        Extra keyword arguments to this function are passed through to the
        corresponding matplotlib method.
        """
        contour_args = self._update_contour_args(contour_args)

        axes = self._check_axes(axes)
        levels = self._process_levels_arg(levels)

        # Pixel indices
        y, x = np.indices(self.data.shape)

        # Prepare a local variable in case we need to mask values
        data = self.data

        # Transform the indices if plotting to a different WCS
        # We do this instead of using the `transform` keyword argument so that Matplotlib does not
        # get confused about the bounds of the contours
        if self.wcs is not axes.wcs:
            if "transform" in contour_args:
                transform_orig = contour_args.pop("transform")
            else:
                transform_orig = axes.get_transform(self.wcs)
            transform = transform_orig - axes.transData  # pixel->pixel transform
            x_1d, y_1d = transform.transform(np.stack([x.ravel(), y.ravel()]).T).T
            x, y = np.reshape(x_1d, x.shape), np.reshape(y_1d, y.shape)

            # Mask out the data array anywhere the coordinate arrays are not finite
            data = np.ma.array(data, mask=~np.logical_and(np.isfinite(x), np.isfinite(y)))

        if fill:
            # Ensure we have more than one level if fill is True
            if len(levels) == 1:
                max_val = np.nanmax(self.data)
                # Ensure the existing level is less than max_val
                if levels[0] < max_val:
                    levels = np.append(levels, max_val)
                else:
                    raise ValueError(
                        f"The provided level ({levels[0]}) is not smaller than the maximum data value ({max_val}). "
                        "Contour level must be smaller than the maximum data value to use `fill=True`.")
            cs = axes.contourf(x, y, data, levels, **contour_args)
        else:
            cs = axes.contour(x, y, data, levels, **contour_args)
        return cs

    def _update_contour_args(self, contour_args):
        """
        Updates ``contour_args`` with values from ``plot_settings``.

        Parameters
        ----------
        contour_args : dict
            A dictionary of arguments to be used for contour plotting.

        Returns
        -------
        dict
            The updated ``contour_args`` dictionary.

        Notes
        -----
        - 'cmap': Set to `None` to avoid the error "ValueError: Either colors or cmap must be None".
        - 'interpolation': Removed because Matplotlib's contour function raises the warning
        "The following kwargs were not used by contour: 'interpolation'".
        - 'origin': If `'origin': 'lower'` is present, it is replaced with `'origin': None`,
        as `None` is the default value for Matplotlib's contour plots.
        """
        plot_settings = self.plot_settings.copy()
        contour_args_copy = contour_args.copy()
        contour_args.update(plot_settings)
        # Define default settings for normal plots and contour-specific updates
        original_plot_defaults = {
            'origin': 'lower',
        }
        default_contour_param = {
            'origin': None,
        }
        # Replace conflicting settings with contour defaults
        for key in original_plot_defaults:
            if key in contour_args and contour_args[key] == original_plot_defaults[key]:
                contour_args[key] = default_contour_param[key]
        # 'cmap' cannot be used for contour plots when levels are not None,
        # which is the case in composite maps.
        contour_args['cmap'] = None
        # custom 'norm' cannot be passed through plot_settings
        contour_args['norm'] = None
        # If 'draw_contour' is used, setting 'norm' and 'cmap' to None ensures the method arguments are applied.
        contour_args.update(contour_args_copy)
        contour_args.pop('interpolation')
        return contour_args


    @peek_show
    def peek(self, draw_limb=False, draw_grid=False,
             colorbar=True, **matplot_args):
        """
        Displays a graphical overview of the data in this object for user evaluation.
        For the creation of plots, users should instead use the `~sunpy.map.GenericMap.plot`
        method and Matplotlib's pyplot framework.

        Parameters
        ----------
        draw_limb : bool
            Whether the solar limb should be plotted.
        draw_grid : bool or `~astropy.units.Quantity`
            Whether solar meridians and parallels are plotted.
            If `~astropy.units.Quantity` then sets degree difference between
            parallels and meridians.
        colorbar : bool
            Whether to display a colorbar next to the plot.
        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting.
        """
        figure = plt.figure()
        axes = wcsaxes_compat.gca_wcs(self.smap.wcs)

        im = self.plot(axes=axes, **matplot_args)

        grid_spacing = None
        # Handle case where draw_grid is actually the grid sapcing
        if isinstance(draw_grid, u.Quantity):
            grid_spacing = draw_grid
            draw_grid = True
        elif not isinstance(draw_grid, bool):
            raise TypeError("draw_grid should be a bool or an astropy Quantity.")

        if colorbar:
            if draw_grid:
                pad = 0.12  # Pad to compensate for ticks and axes labels
            else:
                pad = 0.05  # Default value for vertical colorbar
            colorbar_label = str(self.smap.unit) if self.smap.unit is not None else ""
            figure.colorbar(im, pad=pad).set_label(colorbar_label,
                                                   rotation=0, labelpad=-50, y=-0.02, size=12)

        if draw_limb:
            self.draw_limb(axes=axes)

        if draw_grid:
            if grid_spacing is None:
                self.draw_grid(axes=axes)
            else:
                self.draw_grid(axes=axes, grid_spacing=grid_spacing)

        return figure

    @u.quantity_input
    def plot(self, * , annotate=True, axes=None, title=True, autoalign=False,
             clip_interval: u.percent = None, **imshow_kwargs):
        """
        Plots the map object using matplotlib, in a method equivalent
        to :meth:`~matplotlib.axes.Axes.imshow` using nearest neighbor interpolation.

        Parameters
        ----------
        annotate : `bool`, optional
            If `True`, the data is plotted at its natural scale; with
            title and axis labels.
        axes : `~matplotlib.axes.Axes` or None
            If provided the image will be plotted on the given axes. Else the
            current Matplotlib axes will be used.
        title : `str`, `bool`, optional
            The plot title. If `True`, uses the default title for this map.
        clip_interval : two-element `~astropy.units.Quantity`, optional
            If provided, the data will be clipped to the percentile interval bounded by the two
            numbers.
        autoalign : `bool` or `str`, optional
            If other than `False`, the plotting accounts for any difference between the
            WCS of the map and the WCS of the `~astropy.visualization.wcsaxes.WCSAxes`
            axes (e.g., a difference in rotation angle).  If ``pcolormesh``, this
            method will use :meth:`~matplotlib.axes.Axes.pcolormesh` instead of the
            default :meth:`~matplotlib.axes.Axes.imshow`.  Specifying `True` is
            equivalent to specifying ``pcolormesh``.
        **imshow_kwargs : `dict`
            Any additional imshow arguments are passed to :meth:`~matplotlib.axes.Axes.imshow`.

        Examples
        --------
        >>> # Simple Plot with color bar
        >>> aia.plot()   # doctest: +SKIP
        >>> plt.colorbar()   # doctest: +SKIP
        >>> # Add a limb line and grid
        >>> aia.plot()   # doctest: +SKIP
        >>> aia.draw_limb()   # doctest: +SKIP
        >>> aia.draw_grid()   # doctest: +SKIP

        Notes
        -----
        The ``autoalign`` functionality is computationally intensive.  If the plot will
        be interactive, the alternative approach of preprocessing the map (e.g.,
        de-rotating it) to match the desired axes will result in better performance.

        When combining ``autoalign`` functionality with
        `~sunpy.coordinates.Helioprojective` coordinates, portions of the map that are
        beyond the solar disk may not appear, which may also inhibit Matplotlib's
        autoscaling of the plot limits.  The plot limits can be set manually.
        To preserve the off-disk parts of the map, using the
        :meth:`~sunpy.coordinates.Helioprojective.assume_spherical_screen` context
        manager may be appropriate.
        """
        # Set the default approach to autoalignment
        if autoalign not in [False, True, 'pcolormesh']:
            raise ValueError("The value for `autoalign` must be False, True, or 'pcolormesh'.")
        if autoalign is True:
            autoalign = 'pcolormesh'

        axes = self._check_axes(axes, warn_different_wcs=autoalign is False)

        # Normal plot
        plot_settings = copy.deepcopy(self.plot_settings)
        if 'title' in plot_settings:
            plot_settings_title = plot_settings.pop('title')
        else:
            plot_settings_title = self.smap.latex_name

        # Anything left in plot_settings is given to imshow
        imshow_args = plot_settings
        if annotate:
            if title is True:
                title = plot_settings_title

            if title:
                axes.set_title(title)

            # WCSAxes has unit identifiers on the tick labels, so no need
            # to add unit information to the label
            spatial_units = [None, None]
            ctype = axes.wcs.wcs.ctype

            axes.set_xlabel(axis_labels_from_ctype(ctype[0],
                                                   spatial_units[0]))
            axes.set_ylabel(axis_labels_from_ctype(ctype[1],
                                                   spatial_units[1]))

        # Take a deep copy here so that a norm in imshow_kwargs doesn't get modified
        # by setting it's vmin and vmax
        imshow_args.update(copy.deepcopy(imshow_kwargs))

        if clip_interval is not None:
            if len(clip_interval) == 2:
                clip_percentages = clip_interval.to('%').value
                vmin, vmax = AsymmetricPercentileInterval(*clip_percentages).get_limits(self.smap.data)
            else:
                raise ValueError("Clip percentile interval must be specified as two numbers.")

            imshow_args['vmin'] = vmin
            imshow_args['vmax'] = vmax

        msg = ('Cannot manually specify {0}, as the norm '
               'already has {0} set. To prevent this error set {0} on '
               '`m.plot_settings["norm"]` or the norm passed to `m.plot`.')
        if (norm := imshow_args.get('norm', None)) is not None:
            if 'vmin' in imshow_args:
                if norm.vmin is not None:
                    raise ValueError(msg.format('vmin'))
                norm.vmin = imshow_args.pop('vmin')
            if 'vmax' in imshow_args:
                if norm.vmax is not None:
                    raise ValueError(msg.format('vmax'))
                norm.vmax = imshow_args.pop('vmax')

        if self.smap.mask is None:
            data = self.smap.data
        else:
            data = np.ma.array(np.asarray(self.smap.data), mask=self.smap.mask)

        if autoalign == 'pcolormesh':
            # We have to handle an `aspect` keyword separately
            axes.set_aspect(imshow_args.get('aspect', 1))

            # pcolormesh does not do interpolation
            if imshow_args.get('interpolation', None) not in [None, 'none', 'nearest']:
                warn_user("The interpolation keyword argument is ignored when using autoalign "
                          "functionality.")

            # Remove imshow keyword arguments that are not accepted by pcolormesh
            for item in ['aspect', 'extent', 'interpolation', 'origin']:
                if item in imshow_args:
                    del imshow_args[item]

            imshow_args.setdefault('transform', axes.get_transform(self.smap.wcs))

            # The quadrilaterals of pcolormesh can slightly overlap, which creates the appearance
            # of a grid pattern when alpha is not 1.  These settings minimize the overlap.
            if imshow_args.get('alpha', 1) != 1:
                imshow_args.setdefault('antialiased', True)
                imshow_args.setdefault('linewidth', 0)

            ret = axes.pcolormesh(np.arange(data.shape[1] + 1) - 0.5,
                                  np.arange(data.shape[0] + 1) - 0.5,
                                  data, **imshow_args)
        else:
            ret = axes.imshow(data, **imshow_args)

        wcsaxes_compat.default_wcs_grid(axes)

        # Set current axes/image if pyplot is being used (makes colorbar work)
        for i in plt.get_fignums():
            if axes in plt.figure(i).axes:
                plt.sca(axes)
                plt.sci(ret)

        return ret

    def _check_axes(self, axes, warn_different_wcs=False):
        """
        - If axes is None, get the current Axes object.
        - Error if not a WCSAxes.
        - Return axes.

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            Axes to validate.
        warn_different_wcs : bool
            If `True`, warn if the Axes WCS is different from the Map WCS. This is only used for
            `.plot()`, and can be removed once support is added for plotting a map on a different
            WCSAxes.
        """
        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.smap.wcs)

        if not wcsaxes_compat.is_wcsaxes(axes):
            raise TypeError("The axes need to be an instance of WCSAxes. "
                            "To fix this pass set the `projection` keyword "
                            "to this map when creating the axes.")
        elif warn_different_wcs and not axes.wcs.wcs.compare(self.smap.wcs.wcs, tolerance=0.01):
            warn_user('The map world coordinate system (WCS) is different from the axes WCS. '
                      'The map data axes may not correctly align with the coordinate axes. '
                      'To automatically transform the data to the coordinate axes, specify '
                      '`autoalign=True`.')

        return axes
