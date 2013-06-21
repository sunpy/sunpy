"""A Python MapCube Object"""
from __future__ import absolute_import
#pylint: disable=W0401,W0614,W0201,W0212,W0404

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import copy

from sunpy.map import MapBase
#from .map_factory import Map
#from sunpy.map.sources import *
from sunpy.util import plotting

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

    Attributes
    ----------
    headers : list
        a list of dictionaries containing the original and normalized header tags for the files used to build the MapCube.

    See Also
    --------
    numpy.ndarray Parent class for the MapCube object
    :class:`sunpy.map.Map`
        
    Examples
    --------
    >>> mapcube = sunpy.Map('images/', mapcube=True)
    >>> mapcube[0].show()
    >>> mapcube[3].reference_pixel['x']
    2050.6599120000001
    """
    #pylint: disable=W0613,E1101
    def __init__(self, maps, sortby='date', coalign=False, derotate=False):
        """Creates a new Map instance"""
        
        self._maps = maps

        # Optionally sort data
        if sortby is not None:
            if sortby is 'date':
                self._maps.sort(key=getattr(self, '_sort_by_%s' % sortby)())
            else:
                raise ValueError("Only sort by date is supported")
        
        # Coalignment
        if coalign and hasattr(self, '_coalign_%s' % coalign):
            getattr(self, '_coalign_%s' % coalign)()

        if derotate:
            self._derotate()

    def __getitem__(self, key):
        """Overiding indexing operation"""
        return self._maps[key]
        
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
                    # calculate difference for intersecting pixels
                    # if < min_diff['value'], store new value/offset
                    
            # shift image
            if min_diff['offset'] != (0, 0):
                # shift and clip image

        """
        pass
    
    # Sorting methods
    @classmethod
    def _sort_by_date(cls):
        return lambda m: m.date # maps.sort(key=attrgetter('date'))
    
    def _derotate(self):
        """Derotates the layers in the MapCube"""
        pass
    
    def plot(self, gamma=None, annotate=True, axes=None, controls=True,
             interval=200, resample=False, colorbar=False,
             **ani_args):
        """
        A animation plotting routine that animates each element in the
        MapCube
        
        Parameters
        ----------
        gamma: float
            Gamma value to use for the color map
            
        annotate: bool
            If true, the data is plotted at it's natural scale; with
            title and axis labels.
            
        axes: matplotlib.axes object or None
            If provided the image will be plotted on the given axes. Else the 
            current matplotlib axes will be used.
        
        controls: bool
            Adds play / pause button to the animation
        
        interval: int
            Frame display time in ms.
        
        resample: list or False
            Draws the map at a lower resolution to increase the speed of
            animation. Specify a list as a fraction i.e. [0.25, 0.25] to 
            plot at 1/4 resolution.
        
        colorbar: bool
            Draw a colorbar on the plot.
        
        **ani_args : dict
            Any additional imshow arguments that should be used
            when plotting the image. Passed to 
            sunpy.util.plotting.ControlFuncAnimation
        
        Example
        -------
        cube = MapCube(*maps)
        ani = cube.plot(colorbar=True)        
        plt.show()
        
        #Plot the map at 1/2 original resolution.
        cube = MapCube(*maps)
        ani = cube.plot(resample=[0.5, 0.5], colorbar=True)        
        plt.show()
        """
        
        if not axes:
            axes = plt.gca()
        fig = axes.get_figure()
        
        # Normal plot
        if annotate:
            axes.set_title("%s %s" % (self[0].name, self[0].date))
            
            # x-axis label
            if self[0].coordinate_system['x'] == 'HG':
                xlabel = 'Longitude [%s]' % self[0].units['x']
            else:
                xlabel = 'X-position [%s]' % self[0].units['x']

            # y-axis label
            if self[0].coordinate_system['y'] == 'HG':
                ylabel = 'Latitude [%s]' % self[0].units['y']
            else:
                ylabel = 'Y-position [%s]' % self[0].units['y']
                
            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)

        # Determine extent
        extent = self[0].xrange + self[0].yrange
        
        cmap = copy(self[0].cmap)
        if gamma is not None:
            cmap.set_gamma(gamma)
            
            #make imshow kwargs a dict
        
        kwargs = {'origin':'lower',
                  'cmap':cmap,
                  'norm':self[0].norm(),
                  'extent':extent,
                  'interpolation':'nearest'}
        kwargs.update(ani_args)
        
        im = axes.imshow(self[0], **kwargs)
        
        #Set current image (makes colorbar work)
        plt.sci(im)
        
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.2)
        cbar = plt.colorbar(im,cax)
        
        if resample:
            #This assumes that the maps a homogenous!
            #TODO: Update this!
            resample = np.array(self.shape[1:]) * np.array(resample)
            ani_data = [x.resample(resample) for x in self]
        else:
            ani_data = self
            
        def updatefig(i, *args):
            im = args[0]
            im.set_array(args[2][i])
            im.set_cmap(self[i].cmap)
            im.set_norm(self[i].norm())
            if args[1]:
                axes.set_title("%s %s" % (self[i].name, self[i].date))
        
        ani = plotting.ControlFuncAnimation(fig, updatefig,
                                            frames=xrange(0,self.shape[0]),
                                            fargs=[im,annotate,ani_data],
                                            interval=interval,
                                            blit=False,**ani_args)
        if controls:
            axes, bax1, bax2, bax3 = plotting.add_controls(axes=axes)

            bax1._button.on_clicked(ani._start)
            bax2._button.on_clicked(ani._stop)
            bax3._button.on_clicked(ani._step)
        
        return ani
