"""A Composite Map class

Author: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
from __future__ import absolute_import, print_function, division

import matplotlib.pyplot as plt

import astropy.units as u

from sunpy.map import GenericMap

from sunpy.util import expand_list
from sunpy.extern import six
from sunpy.extern.six.moves import range

__all__ = ['CompositeMap']

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"


class CompositeMap(object):
    """
    CompositeMap(map1 [,map2,..])

    A Composite Map class

    Parameters
    ----------
    args : [`~sunpy.map.Map` | string]
        One or more map of filepaths

    Methods
    -------
    add_map(map, zorder=None, alpha=1, levels=False)
        Adds a map to the CompositeMap
    remove_map(index)
        Removes and returns the map with the given index.
    list_maps()
        Prints a list of the currently included maps.
    get_alpha(index=None)
        Returns the alpha-channel value for a layer in the composite image
    get_zorder(index=None)
        Returns the layering preference (z-order) for a map within the composite.
    get_colors(index=None)
        Returns the colors for a map within the CompositeMap.
    get_norm(index=None)
        Returns the normalization for a map within the CompositeMap.
    get_levels(index=None)
        Returns the list of contour levels for a map within the CompositeMap
    set_norm(self, index, norm)
        Sets the norm for a layer in the composite image.
    set_levels(index, levels, percent=False)
        Sets the contour levels for a layer in the CompositeMap.
    set_colors(index=None, cm)
        Sets the color map for a layer in the CompositeMap.
    set_alpha(index=None, alpha)
        Sets the alpha-channel value for a layer in the CompositeMap.
    set_zorder(index=None, zorder)
        Set the layering preference (z-order) for a map within the CompositeMap.
    plot(figure=None, overlays=None, draw_limb=False, gamma=1.0,
    draw_grid=False, colorbar=True, basic_plot=False,title="SunPy Plot",
    matplot_args)
        Plots the composite map object using matplotlib

    Examples
    --------
    >>> import sunpy.map
    >>> import sunpy.data
    >>> sunpy.data.download_sample_data(overwrite=False)   # doctest: +SKIP
    >>> import sunpy.data.sample
    >>> comp_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE,
    ...                          sunpy.data.sample.EIT_195_IMAGE,
    ...                          composite=True)
    >>> comp_map.add_map(sunpy.map.Map(sunpy.data.sample.RHESSI_IMAGE))
    >>> comp_map.peek()   # doctest: +SKIP

    """
    def __init__(self, *args, **kwargs):
        self._maps = expand_list(args)

        for m in self._maps:
            if not isinstance(m, GenericMap):
                raise ValueError(
                           'CompositeMap expects pre-constructed map objects.')

        # Default alpha and zorder values
        alphas = [1] * len(self._maps)
        zorders = list(range(0, 10 * len(self._maps), 10))
        levels = [False] * len(self._maps)

        # Set z-order and alpha values for the map
        for i, m in enumerate(self._maps):
            m.zorder = zorders[i]
            m.alpha = alphas[i]
            m.levels = levels[i]

    def add_map(self, amap, zorder=None, alpha=1, levels=False):
        """Adds a map to the CompositeMap.

        Parameters
        ----------
        map : `~sunpy.map.GenericMap` or subclass
            Map instance to be added
        zorder : `int`
            The index to use when determining where the map should lie along
            the z-axis; maps with higher z-orders appear above maps with lower
            z-orders.
        alpha : `float`
            Opacity at which the map should be displayed. An alpha value of 0
            results in a fully transparent image while an alpha value of 1
            results in a fully opaque image. Values between result in semi-
            transparent images.
        """
        if zorder is None:
            zorder = max([m.zorder for m in self._maps]) + 10

        amap.zorder = zorder
        amap.alpha = alpha
        amap.levels = levels

        self._maps.append(amap)

    def remove_map(self, index):
        """Removes and returns the map with the given index.

        Parameters
        ----------
        index : `int`
            The index of the map in the composite map.

        Returns
        -------
        `sunpy.map.CompositeMap`
            A composite map with the map indexed by 'index' removed from the
            composite map.
        """
        return self._maps.pop(index)

    def list_maps(self):
        """Prints a list of the currently included maps."""
        print([m.__class__ for m in self._maps])

    def get_map(self, index):
        """Returns the map with given index """
        return self._maps[index]

    def get_alpha(self, index=None):
        """Returns the alpha-channel value for a layer in the composite image"""
        if index is None:
            return [_map.alpha for _map in self._maps]
        else:
            return self._maps[index].alpha

    def get_zorder(self, index=None):
        """Returns the layering preference (z-order) for a map within the
        composite.

        Parameters
        ----------
        index : {`int` | None}
            The index of the map in the composite map.

        Returns
        -------
        {`float` | `list`}
            The layering order (z-order) of the map(s) in the composite
            map.  If None then the layering order of all the maps is returned in
            a list.
        """
        if index is None:
            return [_map.zorder for _map in self._maps]
        else:
            return self._maps[index].zorder

    def get_colors(self, index=None):
        """Returns the colors for a map within the composite map.

        Parameters
        ----------
        index : {`int` | None}
            The index of the map in the composite map.

        Returns
        -------
        {`sunpy.cm` | `list`}
            The colormaps of the map(s) in the composite map.  If None then the
            colormaps of all the maps are returned in a list.
        """

        if index is None:
            return [_map.plot_settings['cmap'] for _map in self._maps]
        else:
            return self._maps[index].plot_settings['cmap']

    def get_mpl_color_normalizer(self, index=None):
        """Returns the color normalizer for a map within the
        composite.

        Parameters
        ----------
        index : {`int` | None}
            The index of the map in the composite map.

        Returns
        -------
        {color normalizer | `list`}
            The color normalizer(s) of the map(s) in the composite map.
            If None then the color normalizers of all the maps are returned in
            a list.
        """
        if index is None:
            return [_map.mpl_color_normalizer for _map in self._maps]
        else:
            return self._maps[index].mpl_color_normalizer

    def get_levels(self, index=None):
        """Returns the list of contour levels for a map within the
        composite.

        Parameters
        ----------
        index : {`int` | None}
            The index of the map in the composite map.

        Returns
        -------
        `list`
            A list of the contour levels of map at index 'index' in the
            composite map.  If index is None, then the contour levels of all
            the maps are returned as a list of lists.
        """
        if index is None:
            return [_map.levels for _map in self._maps]
        else:
            return self._maps[index].levels

    def set_mpl_color_normalizer(self, index, norm):
        """Sets the color normalizer for a layer in the composite image.

        Parameters
        ----------
        index : `int`
            The index of the map in the composite map.

        norm : a color normalizer
            The function used to stretch the color table.

        Returns
        -------
        `~sunpy.map.CompositeMap`
            Sets the color normalizer of the map at index 'index' in the
            composite map to the value given by 'norm'."""
        self._maps[index].mpl_color_normalizer = norm

    def set_levels(self, index, levels, percent=False):
        """
        Sets the contour levels for a layer in the composite image.

        Parameters
        ----------
        index : `int`
            The index of the map in the composite map.

        levels : array-like
            The contour levels.

        percent : `bool`
            If True, the input 'levels' are interpreted as percentages relative
            to the maximum value of the data in layer 'index' of the composite
            map.  If False, the contour levels are set directly from 'levels'.

        Returns
        -------
        `~sunpy.map.CompositeMap`
            A composite map with contour levels 'levels' at layer 'index'.
        """
        if percent is False:
            self._maps[index].levels = levels
        else:
            self._maps[index].levels = [self._maps[index].max()*level/100.0 for level in levels]

    def set_colors(self, index, cm):
        """Sets the color map for a layer in the composite image.

        Parameters
        ----------
        index : `int`
            The index of the map in the composite map.

        cm : a color map
            The contour levels.

        Returns
        -------
        `~sunpy.map.CompositeMap`
            A composite map with colormap 'cm' at layer 'index'.
        """
        self._maps[index].plot_settings['cmap'] = cm

    def set_alpha(self, index, alpha):
        """Sets the alpha-channel value for a layer in the composite image.

        Parameters
        ----------
        index : `int`
            The index of the map in the composite map.

        alpha : `float`
            A float in the range 0 to 1.

        Returns
        -------
        `~sunpy.map.CompositeMap`
            A composite map with alpha-channel value 'alpha' at layer 'index'.
        """
        if 0 <= alpha <= 1:
            self._maps[index].alpha = alpha
        else:
            raise OutOfRangeAlphaValue("Alpha value must be between 0 and 1.")

    def set_zorder(self, index, zorder):
        """Set the layering order (z-order) for a map within the
        composite.

        Parameters
        ----------
        index : `int`
            The index of the map in the composite map.

        zorder : `int`
            The layer order.

        Returns
        -------
        `~sunpy.map.CompositeMap`
            A composite map with the map at layer 'index' having layering order
            'zorder'.
        """
        self._maps[index].zorder = zorder

    def draw_limb(self, index=None, axes=None):
        """Draws a circle representing the solar limb.

        Parameters
        ----------
        index : `int`
            Map index to use to plot limb.

        axes : `matplotlib.axes.Axes` or None
            Axes to plot limb on or None to use current axes.

        Returns
        -------
        `matplotlib.axes.Axes`
        """
        if index is None:
            for i,amap in enumerate(self._maps):
                if hasattr(amap,'rsun_obs'):
                    index = i
                    break

        index_check = hasattr(self._maps[index],'rsun_obs')
        if not index_check or index is None:
            raise ValueError("Specified index does not have all the required attributes to draw limb.")

        return self._maps[index].draw_limb(axes=axes)

    @u.quantity_input(grid_spacing=u.deg)
    def draw_grid(self, index=None, axes=None, grid_spacing=20*u.deg):
        """Draws a grid over the surface of the Sun.

        Parameters
        ----------
        index: int
            Index to determine which map to use to draw grid.

        axes: `~matplotlib.axes.Axes` or None
            Axes to plot limb on or None to use current axes.

        grid_spacing : `float`
            Spacing (in degrees) for longitude and latitude grid.

        Returns
        -------
        `matplotlib.axes.Axes` object
        """
        needed_attrs = ['rsun_meters', 'dsun', 'heliographic_latitude',
                            'heliographic_longitude']
        if index is None:
            for i, amap in enumerate(self._maps):
                if all([hasattr(amap,k) for k in needed_attrs]):
                    index = i
                    break

        index_check = all([hasattr(self._maps[index],k) for k in needed_attrs])
        if not index_check or index is None:
            raise ValueError("Specified index does not have all the required attributes to draw grid.")

        ax = self._maps[index].draw_grid(axes=axes, grid_spacing=grid_spacing)
        return ax

    def plot(self, axes=None, annotate=True, # pylint: disable=W0613
             title="SunPy Composite Plot", **matplot_args):
        """Plots the composite map object using matplotlib

        Parameters
        ----------

        axes: `~matplotlib.axes.Axes` or None
            If provided the image will be plotted on the given axes. Else the
            current matplotlib axes will be used.

        annotate : `bool`
            If true, the data is plotted at it's natural scale; with
            title and axis labels.

        **matplot_args : `dict`
            Matplotlib Any additional imshow arguments that should be used
            when plotting.

        Returns
        -------
        ret : `list`
            List of axes image or quad contour sets that have been plotted.
        """

        # Get current axes
        if not axes:
            axes = plt.gca()

        if annotate:
            # x-axis label
            if self._maps[0].coordinate_system.x == 'HG':
                xlabel = 'Longitude [{lon}]'.format(lon=self._maps[0].spatial_units.x)
            else:
                xlabel = 'X-position [{solx}]'.format(solx=self._maps[0].spatial_units.x)

            # y-axis label
            if self._maps[0].coordinate_system.y == 'HG':
                ylabel = 'Latitude [{lat}]'.format(lat=self._maps[0].spatial_units.y)
            else:
                ylabel = 'Y-position [{soly}]'.format(soly=self._maps[0].spatial_units.y)

            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)

            axes.set_title(title)

        # Define a list of plotted objects
        ret = []
        # Plot layers of composite map
        for m in self._maps:
            # Parameters for plotting
            params = {
                "origin": "lower",
                "extent": list(m.xrange.value) + list(m.yrange.value),
                "cmap": m.plot_settings['cmap'],
                "norm": m.plot_settings['norm'],
                "alpha": m.alpha,
                "zorder": m.zorder,
            }
            params.update(matplot_args)

            if m.levels is False:
                ret.append(axes.imshow(m.data, **params))

            # Use contour for contour data, and imshow otherwise
            if m.levels is not False:
                # Set data with values <= 0 to transparent
                # contour_data = np.ma.masked_array(m, mask=(m <= 0))
                ret.append(axes.contour(m.data, m.levels, **params))
                # Set the label of the first line so a legend can be created
                ret[-1].collections[0].set_label(m.name)

        # Adjust axes extents to include all data
        axes.axis('image')

        # Set current image (makes colorbar work)
        plt.sci(ret[0])
        return ret

    def peek(self, colorbar=True, basic_plot=False, draw_limb=True,
             draw_grid=False, **matplot_args):
        """Displays the map in a new figure.

        Parameters
        ----------
        colorbar : `bool` or `int`
            Whether to display a colorbar next to the plot.
            If specified as an integer a colorbar is plotted for that index.

        basic_plot : `bool`
            If true, the data is plotted by itself at it's natural scale; no
            title, labels, or axes are shown.

        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting.
        """

        # Create a figure and add title and axes
        figure = plt.figure(frameon=not basic_plot)

        # Basic plot
        if basic_plot:
            axes = plt.Axes(figure, [0., 0., 1., 1.])
            axes.set_axis_off()
            figure.add_axes(axes)
            matplot_args.update({'annotate':False})
        else:
            axes = figure.add_subplot(111)

        ret = self.plot(axes=axes,**matplot_args)

        if not isinstance(colorbar, bool) and isinstance(colorbar, int):
            figure.colorbar(ret[colorbar])
        elif colorbar:
            plt.colorbar()
        if draw_limb:
            self.draw_limb(axes=axes)

        if isinstance(draw_grid, bool):
            if draw_grid:
                self.draw_grid(axes=axes)

        elif isinstance(draw_grid, six.integer_types + (float,)):
            self.draw_grid(axes=axes, grid_spacing=draw_grid)
        else:
            raise TypeError("draw_grid should be bool, int, long or float")

        figure.show()


class OutOfRangeAlphaValue(ValueError):
    """Exception to raise when an alpha value outside of the range 0-1 is
    requested.
    """
    pass
