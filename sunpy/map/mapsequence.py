"""A Python MapSequence Object"""
from __future__ import absolute_import, division, print_function
#pylint: disable=W0401,W0614,W0201,W0212,W0404

from copy import deepcopy

import numpy as np
import matplotlib.animation
import numpy.ma as ma

import astropy.units as u

from sunpy.map import GenericMap
from sunpy.visualization.mapsequenceanimator import MapSequenceAnimator
from sunpy.visualization import wcsaxes_compat
from sunpy.visualization import axis_labels_from_ctype
from sunpy.util import expand_list
from sunpy.extern.six.moves import range

__all__ = ['MapSequence']


class MapSequence(object):
    """
    MapSequence

    A series of Maps in a single object.

    Parameters
    ----------
    args : `list`
        A list of Map instances
    sortby : `datetime.datetime`
        Method by which the MapSequence should be sorted along the z-axis.
    derotate : `bool`
        Apply a derotation to the data. Default to False.

    To coalign a mapsequence so that solar features remain on the same pixels,
    please see the "Coalignment of MapSequences" note below.

    Attributes
    ----------
    maps : `list`
        This attribute holds the list of Map instances obtained from parameter args.

    Examples
    --------
    >>> import sunpy.map
    >>> mapsequence = sunpy.map.Map('images/*.fits', sequence=True)   # doctest: +SKIP

    MapSequences can be co-aligned using the routines in sunpy.image.coalignment.
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
                           'MapSequence expects pre-constructed map objects.')

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
        work.  Otherwise, a mapsequence is returned."""

        if isinstance(self.maps[key], GenericMap):
            return self.maps[key]
        else:
            return MapSequence(self.maps[key])

    def __len__(self):
        """Return the number of maps in a mapsequence."""
        return len(self.maps)

    # Sorting methods
    @classmethod
    def _sort_by_date(cls):
        return lambda m: m.date  # maps.sort(key=attrgetter('date'))

    def _derotate(self):
        """Derotates the layers in the MapSequence"""
        pass

    def plot(self, axes=None, resample=None, annotate=True,
             interval=200, plot_function=None, **kwargs):
        """
        A animation plotting routine that animates each element in the
        MapSequence

        Parameters
        ----------
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

        plot_function : function
            A function to be called as each map is plotted. Any variables
            returned from the function will have their ``remove()`` method called
            at the start of the next frame so that they are removed from the plot.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> import matplotlib.animation as animation
        >>> from sunpy.map import Map

        >>> sequence = Map(files, sequence=True)   # doctest: +SKIP
        >>> ani = sequence.plot(colorbar=True)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        Plot the map at 1/2 original resolution

        >>> sequence = Map(files, sequence=True)   # doctest: +SKIP
        >>> ani = sequence.plot(resample=[0.5, 0.5], colorbar=True)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        Save an animation of the MapSequence

        >>> sequence = Map(res, sequence=True)   # doctest: +SKIP

        >>> ani = sequence.plot()   # doctest: +SKIP

        >>> Writer = animation.writers['ffmpeg']   # doctest: +SKIP
        >>> writer = Writer(fps=10, metadata=dict(artist='SunPy'), bitrate=1800)   # doctest: +SKIP

        >>> ani.save('mapsequence_animation.mp4', writer=writer)   # doctest: +SKIP

        Save an animation with the limb at each time step

        >>> def myplot(fig, ax, sunpy_map):
        ...    p = sunpy_map.draw_limb()
        ...    return p
        >>> sequence = Map(files, sequence=True)   # doctest: +SKIP
        >>> ani = sequence.peek(plot_function=myplot)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        """
        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.maps[0].wcs)
        fig = axes.get_figure()

        if not plot_function:
            plot_function = lambda fig, ax, smap: []
        removes = []

        # Normal plot
        def annotate_frame(i):
            axes.set_title("{s.name}".format(s=self[i]))
            axes.set_xlabel(axis_labels_from_ctype(self[i].coordinate_system[0],
                                                   self[i].spatial_units[0]))
            axes.set_ylabel(axis_labels_from_ctype(self[i].coordinate_system[1],
                                                   self[i].spatial_units[1]))

        if resample:
            if self.all_maps_same_shape():
                resample = u.Quantity(self.maps[0].dimensions) * np.array(resample)
                ani_data = [amap.resample(resample) for amap in self.maps]
            else:
                raise ValueError('Maps in mapsequence do not all have the same shape.')
        else:
            ani_data = self.maps

        im = ani_data[0].plot(axes=axes, **kwargs)

        def updatefig(i, im, annotate, ani_data, removes):
            while removes:
                removes.pop(0).remove()

            im.set_array(ani_data[i].data)
            im.set_cmap(ani_data[i].plot_settings['cmap'])

            norm = deepcopy(ani_data[i].plot_settings['norm'])
            # The following explicit call is for bugged versions of Astropy's
            # ImageNormalize
            norm.autoscale_None(ani_data[i].data)
            im.set_norm(norm)

            if wcsaxes_compat.is_wcsaxes(axes):
                im.axes.reset_wcs(ani_data[i].wcs)
                wcsaxes_compat.default_wcs_grid(axes)
            else:
                im.set_extent(np.concatenate((ani_data[i].xrange.value,
                                              ani_data[i].yrange.value)))

            if annotate:
                annotate_frame(i)
            removes += list(plot_function(fig, axes, ani_data[i]))

        ani = matplotlib.animation.FuncAnimation(fig, updatefig,
                                                 frames=list(range(0, len(ani_data))),
                                                 fargs=[im, annotate, ani_data, removes],
                                                 interval=interval,
                                                 blit=False)

        return ani

    def peek(self, resample=None, **kwargs):
        """
        A animation plotting routine that animates each element in the
        MapSequence

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

        plot_function : function
            A function to call to overplot extra items on the map plot.
            For more information see `sunpy.visualization.MapSequenceAnimator`.

        Returns
        -------
        mapsequenceanim : `sunpy.visualization.MapSequenceAnimator`
            A mapsequence animator instance.

        See Also
        --------
        sunpy.visualization.mapsequenceanimator.MapSequenceAnimator

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sunpy.map import Map

        >>> sequence = Map(files, sequence=True)   # doctest: +SKIP
        >>> ani = sequence.peek(colorbar=True)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        Plot the map at 1/2 original resolution

        >>> sequence = Map(files, sequence=True)   # doctest: +SKIP
        >>> ani = sequence.peek(resample=[0.5, 0.5], colorbar=True)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        Plot the map with the limb at each time step

        >>> def myplot(fig, ax, sunpy_map):
        ...    p = sunpy_map.draw_limb()
        ...    return p
        >>> sequence = Map(files, sequence=True)   # doctest: +SKIP
        >>> ani = sequence.peek(plot_function=myplot)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        Decide you want an animation:

        >>> sequence = Map(files, sequence=True)   # doctest: +SKIP
        >>> ani = sequence.peek(resample=[0.5, 0.5], colorbar=True)   # doctest: +SKIP
        >>> mplani = ani.get_animation()   # doctest: +SKIP
        """

        if resample:
            if self.all_maps_same_shape():
                plot_sequence = MapSequence()
                resample = u.Quantity(self.maps[0].dimensions) * np.array(resample)
                for amap in self.maps:
                    plot_sequence.maps.append(amap.resample(resample))
            else:
                raise ValueError('Maps in mapsequence do not all have the same shape.')
        else:
            plot_sequence = self

        return MapSequenceAnimator(plot_sequence, **kwargs)

    def all_maps_same_shape(self):
        """
        Tests if all the maps have the same number pixels in the x and y
        directions.
        """
        return np.all([m.data.shape == self.maps[0].data.shape for m in self.maps])

    def at_least_one_map_has_mask(self):
        """
        Tests if at least one map has a mask.
        """
        return np.any([m.mask is not None for m in self.maps])

    def as_array(self):
        """
        If all the map shapes are the same, their image data is rendered
        into the appropriate numpy object.  If none of the maps have masks,
        then the data is returned as a (ny, nx, nt) ndarray.  If all the maps
        have masks, then the data is returned as a (ny, nx, nt) masked array
        with all the masks copied from each map.  If only some of the maps
        have masked then the data is returned as a (ny, nx, nt) masked array,
        with masks copied from maps as appropriately; maps that do not have a
        mask are supplied with a mask that is full of False entries.
        If all the map shapes are not the same, a ValueError is thrown.
        """
        if self.all_maps_same_shape():
            data = np.swapaxes(np.swapaxes(np.asarray([m.data for m in self.maps]), 0, 1).copy(), 1, 2).copy()
            if self.at_least_one_map_has_mask():
                mask_sequence = np.zeros_like(data, dtype=bool)
                for im, m in enumerate(self.maps):
                    if m.mask is not None:
                        mask_sequence[:, :, im] = m.mask
                return ma.masked_array(data, mask=mask_sequence)
            else:
                return data
        else:
            raise ValueError('Not all maps have the same shape.')

    def all_meta(self):
        """
        Return all the meta objects as a list.
        """
        return [m.meta for m in self.maps]
