# -*- coding: utf-8 -*-

__all__ = ['MapCubeAnimator']

from sunpy.visualization import imageanimator

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
    mapcube: sunpy.map.MapCube
        A MapCube

    annotate: bool
        Annotate the figure with scale and titles

    fig: mpl.figure
        Figure to use

    interval: int
        Animation interval in ms

    colorbar: bool
        Plot colorbar

    Notes
    -----
    Extra keywords are passed to `mapcube[0].plot()` i.e. the `plot()` routine of
    the first map in the cube.
    """
    def __init__(self, mapcube, annotate=True, **kwargs):

        self.mapcube = mapcube
        self.annotate = annotate
        slider_functions = [self.updatefig]
        slider_ranges = [[0,len(mapcube.maps)]]

        imageanimator.BaseFuncAnimator.__init__(self, mapcube.maps, slider_functions,
                                        slider_ranges, **kwargs)

        if annotate:
            self._annotate_plot(0)

    def updatefig(self, val, im, slider):
        i = int(val)
        im.set_array(self.data[i].data)
        im.set_cmap(self.mapcube[i].plot_settings['cmap'])
        im.set_norm(self.mapcube[i].plot_settings['norm'])
        # Having this line in means the plot will resize for non-homogenous
        # maps. However it also means that if you zoom in on the plot bad
        # things happen.
        # im.set_extent(self.mapcube[i].xrange + self.mapcube[i].yrange)
        if self.annotate:
            self._annotate_plot(i)

    def _annotate_plot(self, ind):
        """
        Annotate the image.

        This may overwrite some stuff in `GenericMap.plot()`
        """
        # Normal plot
        self.axes.set_title("{s.name} {s.date!s}".format(s=self.data[ind]))

        # x-axis label
        if self.data[ind].coordinate_system.x == 'HG':
            xlabel = 'Longitude [{lon}]'.format(lon=self.data[ind].units.x)
        else:
            xlabel = 'X-position [{xpos}]'.format(xpos=self.data[ind].units.x)

        # y-axis label
        if self.data[ind].coordinate_system.y == 'HG':
            ylabel = 'Latitude [{lat}]'.format(lat=self.data[ind].units.y)
        else:
            ylabel = 'Y-position [{ypos}]'.format(ypos=self.data[ind].units.y)

        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)

    def plot_start_image(self, ax):
        im = self.mapcube[0].plot(annotate=self.annotate, axes=ax,
                                       **self.imshow_kwargs)
        return im
