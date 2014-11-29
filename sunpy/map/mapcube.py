"""A Python MapCube Object"""
from __future__ import absolute_import
#pylint: disable=W0401,W0614,W0201,W0212,W0404

import numpy as np
import matplotlib.animation
import matplotlib.pyplot as plt

from sunpy.map import GenericMap

from sunpy.visualization.mapcubeanimator import MapCubeAnimator
from sunpy.util import expand_list


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
        """Overiding indexing operation"""
<<<<<<< HEAD
        return self._maps[key]
    
    def coalign(self, method="diff"):
        """ Fine coalign the data"""
        if method == 'diff':
            return self._coalign_diff()
    
    # Coalignment methods
    def _coalign_diff(self):
        """Difference-based coalignment
        
        Coaligns data by minimizing the difference between subsequent images
        before and after shifting the images one to several pixels in each
        direction.
        
        pseudo-code:
        
        for i len(self):
            min_diff = {'value': (), 'offset': (0, 0)} # () is pos infinity
            
            # try shifting 1 pixel in each direction
            for x in (-1, 0, 1):
                for y in (-1, 0, 1):
                    # calculate differenand hasattr(self, '_coalign_%s' % coalign):
            getattr(self, '_coalign_%s' % coalign)()ce for intersecting pixels
                    # if < min_diff['value'], store new value/offset
                    
            # shift image
            if min_diff['offset'] != (0, 0):
                # shift and clip image
=======
        return self.maps[key]

    def __len__(self):
        """Return the number of maps in a mapcube."""
        return len(self.maps)
>>>>>>> 39e7e2bf8d71f8d34816accfbad37f2903702f08

    # Sorting methods
    @classmethod
    def _sort_by_date(cls):
        return lambda m: m.date # maps.sort(key=attrgetter('date'))
<<<<<<< HEAD
    
    def _derotate(self, layer_index=0, clip=True):
        """Derotates the layers in the MapCube relative to a given layer
        
        Inputs
        ------
        layer_index : an integer indicating which layer to calculate the
                      solar derotation offsets against.

        clip : {True, False}
               clip the y,x edges of the each map in the macpube according to
               the maximum solar derotation offsets.
        """
    
        # reference values
        ref_center = self.maps[layer_index].center
        ref_time = self.maps[layer_index].date

        # Number of maps
        nt = len(self.maps)

        # Values of the displacements
        ydiff = np.zeros(nt)
        xdiff = np.zeros(nt)
    
        # Assume all the maps have the same pointing and the Sun rotates in the
        # field of view.  The center of the field of view at the start time gives
        # the initial field of view that then moves at later times.  This initial
        # location has to be rotated forward to get the correct center of the FOV
        # for the next maps.
        for i, m in enumerate(maps):
    
            # Find the center of the field of view
            newx, newy = rot_hpc(ref_center['x'], ref_center['y'],
                                 ref_time, m.date)
    
            # Calculate the displacements in terms of pixels
            if newx is None:
                xdiff[i] = 0.0
            else:
                xdiff[i] = (newx - ref_center['x']) / m.scale['x']
            if newy is None:
                ydiff[i] = 0.0
            else:
                ydiff[i] = (newy - ref_center['y']) / m.scale['y']
    
        # Apply the shifts and mapcube clipping (if requested)
        return _shift_clip(ydiff, xdiff, clip)


    def _shift_clip(self, ydisp, xdisp, clip):
        """
        Helper function that shifts all the maps according to the passed in
        displacements.  Optionally clips the datacube.
        
        Inputs
        ------
        ydisp : array of displacements in the y-direction for the datacube, of
                length nt, where nt is the number of maps in the datacube.
        
        xdisp : array of displacements in the x-direction for the datacube, of
                length nt, where nt is the number of maps in the datacube.
        
        clip : {True, False}
               clip the y,x edges of the datacube according to the maximum
               values in ydisp and xdisp.
        """
        # get the dimensions of the datacube
        # Size of the data
        ny = self.maps[layer_index].shape[0]
        nx = self.maps[layer_index].shape[1]
        nt = len(self.maps)
    
        # Output datacube
        shifted_datacube = np.zeros((ny, nx, nt))
        
        # shift the data cube according to the calculated displacements
        for i, m in enumerate(self.maps):
            shifted_datacube[:, :, i] = shift(m.data, [-ydisp[i], -xdisp[i]])

         # Clip the data if requested
        if clip:
            shifted_datacube = clip_edges(shifted_datacube, ydisp, xdisp)

        # Adjust the mapcube. Adjust the positioning information accordingly.
        for i, m in enumerate(self.maps):
            self.maps[i].meta['xcen'] = self.maps[i].meta['xcen'] + xdisp[i] * m.scale['x']
            self.maps[i].meta['ycen'] = self.maps[i].meta['ycen'] + ydisp[i] * m.scale['y']
            self.maps[i].data = shifted_datacube[:, :, i]
        return self


    def plot(self, gamma=None, annotate=True, axes=None, controls=False,
             interval=200, resample=False, colorbar=False, ani_args={},
             **imshow_args):
=======

    def _derotate(self):
        """Derotates the layers in the MapCube"""
        pass

    def plot(self, gamma=None, axes=None, resample=None, annotate=True,
             interval=200, **kwargs):
>>>>>>> 39e7e2bf8d71f8d34816accfbad37f2903702f08
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
            if self[0].coordinate_system['x'] == 'HG':
                xlabel = 'Longitude [{lon}'.format(lon=self[i].units['x'])
            else:
                xlabel = 'X-position [{xpos}]'.format(xpos=self[i].units['x'])

            # y-axis label
            if self[0].coordinate_system['y'] == 'HG':
                ylabel = 'Latitude [{lat}]'.format(lat=self[i].units['y'])
            else:
                ylabel = 'Y-position [{ypos}]'.format(ypos=self[i].units['y'])

            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)

        if gamma is not None:
            self[0].cmap.set_gamma(gamma)

        if resample:
            #This assumes that the maps are homogenous!
            #TODO: Update this!
            resample = np.array(len(self.maps)-1) * np.array(resample)
            ani_data = [x.resample(resample) for x in self.maps]
        else:
            ani_data = self.maps

        im = ani_data[0].plot(**kwargs)

        def updatefig(i, im, annotate, ani_data):

            im.set_array(ani_data[i].data)
            im.set_cmap(self.maps[i].cmap)
            im.set_norm(self.maps[i].mpl_color_normalizer)
            im.set_extent(self.maps[i].xrange + self.maps[i].yrange)
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
