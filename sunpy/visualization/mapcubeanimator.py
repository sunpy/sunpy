# -*- coding: utf-8 -*-

from copy import deepcopy

from sunpy.visualization import imageanimator, wcsaxes_compat, axis_labels_from_ctype
from sunpy.visualization.wcsaxes_compat import _FORCE_NO_WCSAXES

__all__ = ['MapCubeAnimator']


class MapCubeAnimator(imageanimator.BaseFuncAnimator):
    """
    Create an interactive viewer for a MapCube

    The following keyboard shortcuts are defined in the viewer:

    - 'left': previous step on active slider
    - 'right': next step on active slider
    - 'top': change the active slider up one
    - 'bottom': change the active slider down one
    - 'p': play/pause active slider

    Parameters
    ----------
    mapcube : `sunpy.map.MapCube`
        A MapCube

    annotate : `bool`
        Annotate the figure with scale and titles

    fig : `matplotlib.figure`
        Figure to use

    interval : `int`
        Animation interval in ms

    colorbar : `bool`
        Plot colorbar

    plot_function : function
        A function to call when each map is plotted, the function must have
        the signature `(fig, axes, smap)` where fig and axes are the figure and
        axes objects of the plot and smap is the current frames Map object.
        Any objects returned from this function will have their `remove()` method
        called at the start of the next frame to clear them from the plot.

    Notes
    -----
    Extra keywords are passed to `mapcube[0].plot()` i.e. the `plot()` routine of
    the maps in the cube.
    """

    def __init__(self, mapcube, annotate=True, **kwargs):

        self.mapcube = mapcube
        self.annotate = annotate
        self.user_plot_function = kwargs.pop('plot_function',
                                             lambda fig, ax, smap: [])
        # List of object to remove at the start of each plot step
        self.remove_obj = []
        slider_functions = [self.updatefig]
        slider_ranges = [[0, len(mapcube.maps)]]

        imageanimator.BaseFuncAnimator.__init__(
            self, mapcube.maps, slider_functions, slider_ranges, **kwargs)

        if annotate:
            self._annotate_plot(0)

    def updatefig(self, val, im, slider):
        # Remove all the objects that need to be removed from the
        # plot
        while self.remove_obj:
            self.remove_obj.pop(0).remove()

        i = int(val)
        im.set_array(self.data[i].data)
        im.set_cmap(self.mapcube[i].plot_settings['cmap'])

        norm = deepcopy(self.mapcube[i].plot_settings['norm'])
        # The following explicit call is for bugged versions of Astropy's ImageNormalize
        norm.autoscale_None(self.data[i].data)
        im.set_norm(norm)

        if wcsaxes_compat.is_wcsaxes(im.axes):
            im.axes.reset_wcs(self.mapcube[i].wcs)
            wcsaxes_compat.default_wcs_ticks(im.axes,
                                             self.mapcube[i].spatial_units,
                                             self.mapcube[i].coordinate_system)

        # Having this line in means the plot will resize for non-homogenous
        # maps. However it also means that if you zoom in on the plot bad
        # things happen.
        # im.set_extent(self.mapcube[i].xrange + self.mapcube[i].yrange)
        if self.annotate:
            self._annotate_plot(i)

        self.remove_obj += list(
            self.user_plot_function(self.fig, self.axes, self.mapcube[i]))

    def _annotate_plot(self, ind):
        """
        Annotate the image.

        This may overwrite some stuff in `GenericMap.plot()`
        """
        # Normal plot
        self.axes.set_title("{s.name}".format(s=self.data[ind]))

        self.axes.set_xlabel(axis_labels_from_ctype(self.data[ind].coordinate_system[0],
                                                    self.data[ind].spatial_units[0]))
        self.axes.set_ylabel(axis_labels_from_ctype(self.data[ind].coordinate_system[1],
                                                    self.data[ind].spatial_units[1]))

    def _get_main_axes(self):
        """
        Create an axes which is wcsaxes if we have that...
        """
        if not _FORCE_NO_WCSAXES:
            return self.fig.add_subplot(111, projection=self.mapcube[0].wcs)
        else:
            return self.fig.add_subplot(111)

    def plot_start_image(self, ax):
        im = self.mapcube[0].plot(
            annotate=self.annotate, axes=ax, **self.imshow_kwargs)
        self.remove_obj += list(
            self.user_plot_function(self.fig, self.axes, self.mapcube[0]))
        return im
