"""A Python MapCube Object"""
from __future__ import absolute_import
#pylint: disable=W0401,W0614,W0201,W0212,W0404

import numpy as np
import matplotlib.animation
import matplotlib.pyplot as plt

from sunpy.map import GenericMap

from sunpy.visualization.mapcubeanimator import MapCubeAnimator
from sunpy.util import expand_list

# Mapcube co-alignment functions
from sunpy.image.coalignment import default_fmap_function, calculate_shift, clip_edges
from scipy.ndimage.interpolation import shift

__all__ = ['MapCube']

class MapCube(object):
    """
    MapCube(input)

    A series of spatially aligned Maps.

    Parameters
    ----------
    args : {List}
        A list of Map instances
    sortby : {"date", None}
        Method by which the MapCube should be sorted along the z-axis.
    derotate : {None}
        Apply a derotation to the data (Not Implemented)
    coalign : {None}
        Apply fine coalignment to the data (NOTE: requires the installation of
        scikit-image.)

    Attributes
    ----------
    maps : {List}
        This attribute holds the list of Map instances obtained from parameter args.

    Examples
    --------
    >>> mapcube = sunpy.Map('images/', mapcube=True)
    >>> mapcube[0].plot()
    >>> mapcube[3].reference_pixel['x']
    2050.6599120000001
    """
    #pylint: disable=W0613,E1101
    def __init__(self, *args, **kwargs):
        """Creates a new Map instance"""

        # Hack to get around Python 2.x not backporting PEP 3102.
        sortby = kwargs.pop('sortby', 'date')
        coalign = kwargs.pop('coalign', False)
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

        # Coalignment
        if coalign:
            if coalign == 'diff':
                self.coalign("diff")
            else:
                raise ValueError("That coalignment method is not supported")

        if derotate:
            self._derotate()

    def __getitem__(self, key):
        """Overiding indexing operation"""
        return self.maps[key]

    def coalign(self, method="match_template", **kwargs):
        """ Fine coalign the data"""
        if method == 'match_template':
            return self._coalign_by_match_template(**kwargs)
        else:
            raise ValueError("Only 'match_template' coalignment method is supported at present.")

    # Coalignment by matching a template
    def _coalign_by_match_template(self, layer_index=0, func=default_fmap_function,
                                  clip=True, template=None):        
        """
        Co-register the layers in a mapcube according to a template taken from
        that mapcube.  This method requires that scikit-image be installed.
    
        Input
        -----
        self : a mapcube of shape (ny, nx, nt), where nt is the number of
             layers in the mapcube.
    
        layer_index : the layer in the mapcube from which the template will be
                      extracted.
    
        func: a function which is applied to the data values before the
              coalignment method is applied.  This can be useful in coalignment,
              because it is sometimes better to co-align on a function of the data
              rather than the data itself.  The calculated shifts are applied to
              the original data.  Useful functions to consider are the log of the
              image data, or 1 / data. The function is of the form func = F(data).
              The default function ensures that the data are floats.
    
        clip : clip off x, y edges in the datacube that are potentially affected
                by edges effects.
        
        template: {None, Map, ndarray}
                  The template used in the matching.  The template can be
                  another SunPy map, or a numpy ndarray.
    
        Output
        ------
        datacube : the input datacube each layer having been co-registered against
                   the template.
    
        """
        # Size of the data
        ny = self.maps[layer_index].shape[0]
        nx = self.maps[layer_index].shape[1]
        nt = len(self.maps)
    
        # Storage for the shifted data and the pixel shifts
        shifted_datacube = np.zeros((ny, nx, nt))
        xshift_keep = np.zeros((nt))
        yshift_keep = np.zeros((nt))
    
        # Calculate a template.  If no template is passed then define one
        # from the the index layer.
        if template is None:
            tplate = self.maps[layer_index].data[ny / 4: 3 * ny / 4,
                                             nx / 4: 3 * nx / 4]
        elif isinstance(template, GenericMap):
            tplate = template.data
        elif isinstance(template, np.ndarray):
            tplate = template
        else:
            raise ValueError('Invalid template.')
        
        # Apply the function to the template
        tplate = func(tplate)

        # Match the template and calculate shifts
        for i, m in enumerate(self.maps):
            # Get the next 2-d data array
            this_layer = func(m.data)
    
            # Calculate the y and x shifts in pixels
            yshift, xshift = calculate_shift(this_layer, tplate)
    
            # Keep shifts in pixels
            yshift_keep[i] = yshift
            xshift_keep[i] = xshift
    
        # Calculate shifts relative to the template layer
        yshift_keep = yshift_keep - yshift_keep[layer_index]
        xshift_keep = xshift_keep - xshift_keep[layer_index]
    
        # Shift the data
        for i, m in enumerate(self.maps):
            shifted_datacube[:, :, i] = shift(m.data, [-yshift_keep[i], -xshift_keep[i]])
    
        # Clip the data if requested
        if clip:
            shifted_datacube = clip_edges(shifted_datacube, yshift_keep, xshift_keep)
    
        # Create a new mapcube.  Adjust the positioning information accordingly.
        for i, m in enumerate(self.maps):
            self.maps[i].meta['xcen'] = self.maps[i].meta['xcen']  + xshift_keep[i] * m.scale['x']
            self.maps[i].meta['ycen'] = self.maps[i].meta['ycen']  + yshift_keep[i] * m.scale['y']
            self.maps[i].data = shifted_datacube[:, :, i]
        return self

    # Sorting methods
    @classmethod
    def _sort_by_date(cls):
        return lambda m: m.date # maps.sort(key=attrgetter('date'))

    def _derotate(self):
        """Derotates the layers in the MapCube"""
        pass

    def plot(self, gamma=None, axes=None, resample=None, annotate=True,
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
            axes.set_title("%s %s" % (self[i].name, self[i].date))

            # x-axis label
            if self[0].coordinate_system['x'] == 'HG':
                xlabel = 'Longitude [%s]' % self[i].units['x']
            else:
                xlabel = 'X-position [%s]' % self[i].units['x']

            # y-axis label
            if self[0].coordinate_system['y'] == 'HG':
                ylabel = 'Latitude [%s]' % self[i].units['y']
            else:
                ylabel = 'Y-position [%s]' % self[i].units['y']

            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)

        if gamma is not None:
            self[0].cmap.set_gamma(gamma)

        if resample:
            #This assumes that the maps a homogenous!
            #TODO: Update this!
            resample = np.array(len(self.maps)-1) * np.array(resample)
            ani_data = [x.resample(resample) for x in self]
        else:
            ani_data = self

        im = ani_data[0].plot(**kwargs)

        def updatefig(i, im, annotate, ani_data):

            im.set_array(ani_data[i].data)
            im.set_cmap(self[i].cmap)
            im.set_mpl_color_normalizer(self[i].mpl_color_normalizer)
            im.set_extent(self.xrange + self.yrange)
            if annotate:
                annotate_frame(i)

        ani = matplotlib.animation.FuncAnimation(fig, updatefig,
                                            frames=range(0,len(self.maps)),
                                            fargs=[im,annotate,ani_data],
                                            interval=interval,
                                            blit=False)

        return ani

    def peek(self, gamma=None, resample=None, **kwargs):
        """
        A animation plotting routine that animates each element in the
        MapCube

        Parameters
        ----------
        gamma: float
            Gamma value to use for the color map

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

        if gamma is not None:
            self[0].cmap.set_gamma(gamma)

        if resample:
            #This assumes that the maps a homogenous!
            #TODO: Update this!
            resample = np.array(len(self.maps)-1) * np.array(resample)
            for amap in self.maps:
                amap.resample(resample)

        return MapCubeAnimator(self, **kwargs)
