"""
This module provides a way to animate `~sunpy.map.MapSequence`.
"""
from copy import deepcopy

from mpl_animators import BaseFuncAnimator

from sunpy.visualization import axis_labels_from_ctype, wcsaxes_compat

__all__ = ['MapSequenceAnimator']


class MapSequenceAnimator(BaseFuncAnimator):
    """
    Create an interactive viewer for a `~sunpy.map.MapSequence`.

    The following keyboard shortcuts are defined in the viewer:

    * 'left': previous step on active slider.
    * 'right': next step on active slider.
    * 'top': change the active slider up one.
    * 'bottom': change the active slider down one.
    * 'p': play/pause active slider.

    Parameters
    ----------
    mapsequence : `sunpy.map.MapSequence`
        A `~sunpy.map.MapSequence`.
    annotate : `bool`
        Annotate the figure with scale and titles.
    fig : `matplotlib.figure.Figure`
        Figure to use.
    interval : `int`
        Animation interval in milliseconds.
    colorbar : `bool`
        Plot colorbar.
    plot_function : `function`
        A function to call when each `~sunpy.map.Map` is plotted, the function must have
        the signature ``(fig, axes, smap)`` where ``fig`` and ``axes`` are the figure and
        axes objects of the plot and ``smap`` is the current frame's `~sunpy.map.Map` object.
        Any objects returned from this function will have their ``remove()`` method
        called at the start of the next frame to clear them from the plot.

    Notes
    -----
    Extra keywords are passed to ``mapsequence[0].plot()`` i.e. the ``plot()`` routine of
    the maps in the sequence.
    """

    def __init__(self, mapsequence, annotate=True, **kwargs):

        self.mapsequence = mapsequence
        self.annotate = annotate
        self.user_plot_function = kwargs.pop('plot_function',
                                             lambda fig, ax, smap: [])
        # List of object to remove at the start of each plot step
        self.remove_obj = []
        slider_functions = [self.updatefig]
        slider_ranges = [[0, len(mapsequence.maps)]]

        super().__init__(mapsequence.maps, slider_functions, slider_ranges, **kwargs)

        if annotate:
            self._annotate_plot(0)

    def updatefig(self, val, im, slider):
        # Remove all the objects that need to be removed from the
        # plot
        while self.remove_obj:
            self.remove_obj.pop(0).remove()
        i = int(val)
        im.set_array(self.data[i].data)
        im.set_cmap(self.mapsequence[i].plot_settings.get('cmap', "grey"))
        if norm := self.mapsequence[i].plot_settings.get('norm'):
            im.set_norm(deepcopy(norm))
        if wcsaxes_compat.is_wcsaxes(im.axes):
            im.axes.reset_wcs(self.mapsequence[i].wcs)
        # Having this line in means the plot will resize for non-homogeneous
        # maps. However it also means that if you zoom in on the plot bad
        # things happen.
        # im.set_extent(self.mapsequence[i].xrange + self.mapsequence[i].yrange)
        if self.annotate:
            self._annotate_plot(i)
        self.remove_obj += list(
            self.user_plot_function(self.fig, self.axes, self.mapsequence[i]))

    def _annotate_plot(self, ind):
        """
        Annotate the image.

        This may overwrite some stuff in `sunpy.map.GenericMap.plot`
        """
        # Normal plot
        self.axes.set_title(f"{self.data[ind].name}")

        self.axes.set_xlabel(axis_labels_from_ctype(self.data[ind].coordinate_system[0],
                                                    self.data[ind].spatial_units[0]))
        self.axes.set_ylabel(axis_labels_from_ctype(self.data[ind].coordinate_system[1],
                                                    self.data[ind].spatial_units[1]))

    def _setup_main_axes(self):
        """
        Create an axes which is a `~astropy.visualization.wcsaxes.WCSAxes`.
        """
        if self.axes is None:
            self.axes = self.fig.add_subplot(111, projection=self.mapsequence[0].wcs)

    def plot_start_image(self, ax):
        im = self.mapsequence[0].plot(
            annotate=self.annotate, axes=ax, **self.imshow_kwargs)
        self.remove_obj += list(
            self.user_plot_function(self.fig, self.axes, self.mapsequence[0]))
        return im
