"""A Python MapSequence Object"""

import html
import textwrap
import warnings
import webbrowser
from copy import deepcopy
from tempfile import NamedTemporaryFile

import matplotlib.animation
import numpy as np
import numpy.ma as ma

import astropy.units as u

from sunpy.map import GenericMap
from sunpy.util import SunpyUserWarning, expand_list
from sunpy.visualization import axis_labels_from_ctype, wcsaxes_compat
from sunpy.visualization.animator.mapsequenceanimator import MapSequenceAnimator

__all__ = ['MapSequence']


class MapSequence:
    """
    MapSequence

    A series of Maps in a single object.

    Parameters
    ----------
    args : `list`
        A list of Map instances
    sortby : { "date" | `None`}
        Method by which the MapSequence should be sorted along the z-axis.
        Defaults to sorting by: "date" and is the only supported sorting strategy.
        Passing `None` will disable sorting.
    derotate : `bool`
        Apply a derotation to the data. Default to False.

    Attributes
    ----------
    maps : `list`
        This attribute holds the list of Map instances obtained from parameter args.

    Notes
    -----
    To coalign a mapsequence so that solar features remain on the same pixels,
    please see the "Coalignment of MapSequences" note below.

    Examples
    --------
    >>> import sunpy.map
    >>> mapsequence = sunpy.map.Map('images/*.fits', sequence=True)   # doctest: +SKIP

    MapSequences can be co-aligned using the routines in sunpy.image.coalignment.
    """

    def __init__(self, *args, sortby='date', derotate=False, **kwargs):
        """Creates a new Map instance"""

        self.maps = expand_list(args)

        for m in self.maps:
            if not isinstance(m, GenericMap):
                raise ValueError('MapSequence expects pre-constructed map objects.')

        # Optionally sort data
        if sortby is not None:
            if sortby == 'date':
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

    def __repr__(self):
        names = set([m.__class__.__name__ for m in self.maps])
        return (object.__repr__(self) + "\n" +
                f'MapSequence of {len(self.maps)} elements, with maps from {", ".join(names)}')

    def _repr_html_(self):
        nmaps = len(self)

        # Output a warning about rendering time if there are more than 9 Maps
        if nmaps > 9:
            warnings.warn(f"Rendering the summary for a MapSequence of {nmaps} Maps "
                          "may take a while.", SunpyUserWarning)

        # Assemble the individual HTML repr from each Map, all hidden initally
        repr_list = [f"<div style='display: none' index={i}>{m._repr_html_()}</div>"
                     for i, m in enumerate(self.maps)]

        # Unhide the first Map
        repr_list_html = "\n".join(repr_list).replace('display: none', 'display: ', 1)

        # Return HTML with Javascript-powered buttons
        # To avoid potential conflicts, the Javascript code does not use any user-defined functions
        return textwrap.dedent(f"""\
            <pre>{html.escape(self.__repr__())}</pre>
            <form cur_index=0 max_index={nmaps - 1}>
                <!-- Button to decrement index (always starts disabled) -->
                <input type=button value='&larr;' style='font-weight: bold' disabled onClick='
                    var form = this.parentElement;

                    // Decrement index if allowed
                    var cur_index = Math.max(
                                        parseInt(form.getAttribute("cur_index")) - 1,
                                        0
                                    );
                    form.setAttribute("cur_index", cur_index);

                    // Enable the decrement button if and only if this is not the first Map
                    form.children[0].disabled = (cur_index == 0);

                    // Always enable the increment button (because we just decremented)
                    form.children[1].disabled = false;

                    // Update string (which is children[2] of the form)
                    form.children[2].innerHTML = "Map at index " + cur_index;

                    // Update visibilities to show only the current index
                    // This avoids for...of syntax to retain support for ES5 browsers (e.g., IE11)
                    var array = Array.prototype.slice.call(form.lastElementChild.children);
                    array.forEach(function (elem)
                        {{
                            var form = elem.parentElement.parentElement;
                            elem.style.display = (elem.getAttribute("index") ==
                                                      form.getAttribute("cur_index") ? "" : "none"
                                                 );
                        }}
                    );
                '/>

                <!-- Button to increment index (starts enabled if there is more than one Map) -->
                <input type=button value='&rarr;' style='font-weight: bold'
                    {"" if nmaps > 1 else "disabled"} onClick='

                    var form = this.parentElement;

                    // Increment index if allowed
                    var cur_index = Math.min(
                                        parseInt(form.getAttribute("cur_index")) + 1,
                                        form.getAttribute("max_index")
                                    );
                    form.setAttribute("cur_index", cur_index);

                    // Always enable the decrement button (because we just incremented)
                    form.children[0].disabled = false;

                    // Enable the increment button if and only if this is not the last Map
                    form.children[1].disabled = (cur_index == form.getAttribute("max_index"));

                    // Update string (which is children[2] of the form)
                    form.children[2].innerHTML = "Map at index " + cur_index;

                    // Update visibilities to show only the current index
                    // This avoids for...of syntax to retain support for ES5 browsers (e.g., IE11)
                    var array = Array.prototype.slice.call(form.lastElementChild.children);
                    array.forEach(function (elem)
                        {{
                            var form = elem.parentElement.parentElement;
                            elem.style.display = (elem.getAttribute("index") ==
                                                      form.getAttribute("cur_index") ? "" : "none"
                                                 );
                        }}
                    );

                '/>

                <!-- This string is updated as the index is changed -->
                <span>Map at index 0</span>

                <!-- This element is at the end so that lastElementChild will point to it -->
                <div>
                    {repr_list_html}
                </div>
            </form>""")

    def quicklook(self):
        """
        Display a quicklook summary of the MapSequence instance using the default web browser.

        Click on the |larr| and |rarr| buttons to step through the individual maps.

        .. |larr|   unicode:: U+02190 .. LEFTWARDS ARROW
        .. |rarr|   unicode:: U+02192 .. RIGHTWARDS ARROW

        Notes
        -----
        The image colormap uses
        `histogram equalization <https://en.wikipedia.org/wiki/Histogram_equalization>`__.

        Interactive elements require Javascript support to be enabled in the web browser.

        Examples
        --------
        >>> from sunpy.map import Map
        >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
        >>> seq = Map(sunpy.data.sample.HMI_LOS_IMAGE,
        ...           sunpy.data.sample.AIA_1600_IMAGE,
        ...           sunpy.data.sample.EIT_195_IMAGE,
        ...           sequence=True)  # doctest: +REMOTE_DATA
        >>> seq.quicklook()  # doctest: +SKIP

        (which will open the following content in the default web browser)

        .. generate:: html
            :html_border:

            from sunpy.map import Map
            import sunpy.data.sample
            seq = Map(sunpy.data.sample.HMI_LOS_IMAGE,
                      sunpy.data.sample.AIA_1600_IMAGE,
                      sunpy.data.sample.EIT_195_IMAGE,
                      sequence=True)
            print(seq._repr_html_())

        """
        with NamedTemporaryFile('w', delete=False, prefix='sunpy.map.', suffix='.html') as f:
            url = 'file://' + f.name
            f.write(textwrap.dedent(f"""\
                <html>
                    <title>Quicklook summary for {html.escape(object.__repr__(self))}</title>
                    <body>{self._repr_html_()}</body>
                </html>"""))
        webbrowser.open_new_tab(url)

    # Sorting methods
    @classmethod
    def _sort_by_date(cls):
        return lambda m: m.date  # maps.sort(key=attrgetter('date'))

    def _derotate(self):
        """Derotates the layers in the MapSequence"""
        raise NotImplementedError("This functionality has not yet been implemented.")

    def plot(self, axes=None, resample=None, annotate=True,
             interval=200, plot_function=None, **kwargs):
        """
        A animation plotting routine that animates each element in the
        MapSequence

        Parameters
        ----------
        axes: matplotlib.axes.Axes
            axes to plot the animation on, if none uses current axes

        resample: list
            Draws the map at a lower resolution to increase the speed of
            animation. Specify a list as a fraction i.e. [0.25, 0.25] to
            plot at 1/4 resolution.
            [Note: this will only work where the map arrays are the same size]

        annotate: bool
            Annotate the figure with scale and titles

        interval: int
            Animation interval in ms

        plot_function : function
            A function to be called as each map is plotted.
            For more information see `sunpy.visualization.animator.MapSequenceAnimator`.

        Returns
        -------
        `matplotlib.animation.FuncAnimation`
            A FuncAnimation instance.

        See Also
        --------
        `sunpy.visualization.animator.MapSequenceAnimator`

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
            def plot_function(fig, ax, smap): return []
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
                bl = ani_data[i]._get_lon_lat(ani_data[i].bottom_left_coord)
                tr = ani_data[i]._get_lon_lat(ani_data[i].top_right_coord)
                x_range = list(u.Quantity([bl[0], tr[0]]).to(ani_data[i].spatial_units[0]).value)
                y_range = list(u.Quantity([bl[1], tr[1]]).to(ani_data[i].spatial_units[1]).value)

                im.set_extent(np.concatenate((x_range.value, y_range.value)))

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
        fig: matplotlib.figure.Figure
            Figure to use to create the explorer

        resample: list
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
            For more information see `sunpy.visualization.animator.MapSequenceAnimator`.

        Returns
        -------
        mapsequenceanim : `sunpy.visualization.animator.MapSequenceAnimator`

        See Also
        --------
        sunpy.visualization.animator.MapSequenceAnimator

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
            data = np.swapaxes(np.swapaxes(np.asarray(
                [m.data for m in self.maps]), 0, 1).copy(), 1, 2).copy()
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
