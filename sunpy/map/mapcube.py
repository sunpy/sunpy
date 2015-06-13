"""A Python MapCube Object"""
from __future__ import absolute_import
#pylint: disable=W0401,W0614,W0201,W0212,W0404

import numpy as np
import matplotlib.animation
from sunpy.visualization import wcsaxes_compat
import matplotlib.pyplot as plt

from sunpy.map import GenericMap

from sunpy.visualization.mapcubeanimator import MapCubeAnimator
from sunpy.util import expand_list

__all__ = ['MapCube']

class MapCube(object):
    """
    MapCube

    A series of spatially aligned Maps.

    Parameters
    ----------
    args : {List}
        A list of Map instances
    sortby : {"date", None}
        Method by which the MapCube should be sorted along the z-axis.
    derotate : {None}
        Apply a derotation to the data (Not Implemented)

    To coalign a mapcube so that solar features remain on the same pixels,
    please see the "Coalignment of mapcubes" note below.

    Attributes
    ----------
    maps : {List}
        This attribute holds the list of Map instances obtained from parameter args.

    Examples
    --------
    >>> mapcube = sunpy.map.Map('images/*.fits', cube=True)

    Mapcubes can be co-aligned using the routines in sunpy.image.coalignment.
    """
    #pylint: disable=W0613,E1101
    def __init__(self, *args, **kwargs):
        """Creates a new Map instance"""

        # Hack to get around Python 2.x not backporting PEP 3102.
        sortby = kwargs.pop('sortby', 'date')
        derotate = kwargs.pop('derotate', False)

        self.maps = expand_list(args)

        for m in self.maps:
            if not isinstance(m, GenericMap):
                raise ValueError(
                           'CompositeMap expects pre-constructed map objects.')

        # Optionally sort data
        if sortby is not None:
            if sortby is 'date':
                self.maps.sort(key=self._sort_by_date())
            else:
                raise ValueError("Only sort by date is supported")

        if derotate:
            self._derotate()

    def __getitem__(self, key):
        """Overriding indexing operation.  If the key results in a single map,
        then a map object is returned.  This allows functions like enumerate to
        work.  Otherwise, a mapcube is returned."""

        if isinstance(self.maps[key], GenericMap):
            return self.maps[key]
        else:
            return MapCube(self.maps[key])

    def __len__(self):
        """Return the number of maps in a mapcube."""
        return len(self.maps)

    # Sorting methods
    @classmethod
    def _sort_by_date(cls):
        return lambda m: m.date # maps.sort(key=attrgetter('date'))

    def _derotate(self):
        """Derotates the layers in the MapCube"""
        pass

    def plot(self, axes=None, resample=None, annotate=True,
             interval=200, **kwargs):
        """
        A animation plotting routine that animates each element in the
        MapCube

        Parameters
        ----------
        gamma: float
            Gamma value to use for the color map

        axes: mpl axes
            axes to plot the animation on, if none uses current axes

        resample: list or False
            Draws the map at a lower resolution to increase the speed of
            animation. Specify a list as a fraction i.e. [0.25, 0.25] to
            plot at 1/4 resolution.
            [Note: this will only work where the map arrays are the same size]

        annotate: bool
            Annotate the figure with scale and titles

        interval: int
            Animation interval in ms

        Examples
        --------
        >>> cube = sunpy.Map(files, cube=True)
        >>> ani = cube.plot(colorbar=True)
        >>> plt.show()

        Plot the map at 1/2 original resolution

        >>> cube = sunpy.Map(files, cube=True)
        >>> ani = cube.plot(resample=[0.5, 0.5], colorbar=True)
        >>> plt.show()

        Save an animation of the MapCube

        >>> cube = sunpy.Map(res, cube=True)

        >>> ani = cube.plot(controls=False)

        >>> Writer = animation.writers['ffmpeg']
        >>> writer = Writer(fps=10, metadata=dict(artist='SunPy'), bitrate=1800)

        >>> ani.save('mapcube_animation.mp4', writer=writer)
        """
        if not axes:
            axes = plt.gca()
        fig = axes.get_figure()

        # Normal plot
        def annotate_frame(i):
            axes.set_title("{s.name} {s.date!s}".format(s=self[i]))

            # x-axis label
            if self[0].coordinate_system.x == 'HG':
                xlabel = 'Longitude [{lon}'.format(lon=self[i].units.x)
            else:
                xlabel = 'X-position [{xpos}]'.format(xpos=self[i].units.x)

            # y-axis label
            if self[0].coordinate_system.y == 'HG':
                ylabel = 'Latitude [{lat}]'.format(lat=self[i].units.y)
            else:
                ylabel = 'Y-position [{ypos}]'.format(ypos=self[i].units.y)

            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)

        if resample:
            #This assumes that the maps are homogenous!
            #TODO: Update this!
            resample = np.array(len(self.maps)-1) * np.array(resample)
            ani_data = [x.resample(resample) for x in self.maps]
        else:
            ani_data = self.maps

        im = ani_data[0].plot(axes=axes, **kwargs)

        def updatefig(i, im, annotate, ani_data):

            im.set_array(ani_data[i].data)
            im.set_cmap(self.maps[i].plot_settings['cmap'])
            im.set_norm(self.maps[i].plot_settings['norm'])
            im.set_extent(self.maps[i].xrange + self.maps[i].yrange)
            if annotate:
                annotate_frame(i)

        ani = matplotlib.animation.FuncAnimation(fig, updatefig,
                                                frames=range(0,len(self.maps)),
                                                fargs=[im,annotate,ani_data],
                                                interval=interval,
                                                blit=False)

        return ani

    def peek(self, resample=None, **kwargs):
        """
        A animation plotting routine that animates each element in the
        MapCube

        Parameters
        ----------
        fig: mpl.figure
            Figure to use to create the explorer

        resample: list or False
            Draws the map at a lower resolution to increase the speed of
            animation. Specify a list as a fraction i.e. [0.25, 0.25] to
            plot at 1/4 resolution.
            [Note: this will only work where the map arrays are the same size]

        annotate: bool
            Annotate the figure with scale and titles

        interval: int
            Animation interval in ms

        colorbar: bool
            Plot colorbar

        Returns
        -------
        Returns a MapCubeAnimator object

        See Also
        --------
        sunpy.visualization.mapcubeanimator.MapCubeAnimator

        Examples
        --------
        >>> cube = sunpy.Map(files, cube=True)
        >>> ani = cube.plot(colorbar=True)
        >>> plt.show()

        Plot the map at 1/2 original resolution

        >>> cube = sunpy.Map(files, cube=True)
        >>> ani = cube.plot(resample=[0.5, 0.5], colorbar=True)
        >>> plt.show()

        Decide you want an animation:

        >>> cube = sunpy.Map(files, cube=True)
        >>> ani = cube.plot(resample=[0.5, 0.5], colorbar=True)
        >>> mplani = ani.get_animation()
        """

        if resample:
            if self.all_maps_same_shape():
                resample = np.array(len(self.maps) - 1) * np.array(resample)
                for amap in self.maps:
                    amap.resample(resample)
            else:
                raise ValueError('Maps in mapcube do not all have the same shape.')

        return MapCubeAnimator(self, **kwargs)

    def all_maps_same_shape(self):
        """
        Tests if all the maps have the same number pixels in the x and y
        directions.
        """
        return np.all([m.data.shape == self.maps[0].data.shape for m in self.maps])

    def as_array(self):
        """
        If all the map shapes are the same, their image data is copied
        into a single single ndarray. The ndarray is ordered as (ny, nx, nt).
        Otherwise, a ValueError is thrown.
        """
        if self.all_maps_same_shape():
            return np.swapaxes(np.swapaxes(np.asarray([m.data for m in self.maps]), 0, 1).copy(), 1, 2).copy()
        else:
            raise ValueError('Not all maps have the same shape.')

    def all_meta(self):
        """
        Return all the meta objects as a list.
        """
        return [m.meta for m in self.maps]
