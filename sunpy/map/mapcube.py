"""A Python MapCube Object"""
from __future__ import absolute_import
#pylint: disable=W0401,W0614,W0201,W0212,W0404

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
import matplotlib.pyplot as plt

from sunpy.map import Map
from sunpy.map.sources import *
import numpy as np

__all__ = ['MapCube']

# (https://github.com/sunpy/sunpy/issues/397)
# 2011/04/13: Should Map be broken up into Map and MapHeader classes? This way
# mapped header values can be used in MapCube without having to keep extra
# copies of the data..
#
class MapCube(np.ndarray):
    """
    MapCube(input)
    
    A spatially-aware data array based on the SolarSoft Map object.
    Reads in the files at the specified location, stores their headers, and
    creates a 3d array from their contents.

    Parameters
    ----------
    args : {string | Map}* 
        Map instances or filepaths from which MapCube should be built.
    sortby : {"date"}
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
    >>> mapcube = sunpy.make_map('images/')
    >>> mapcube[0].show()
    >>> mapcube[3].reference_pixel['x']
    2050.6599120000001
    """
    def __new__(cls, *args, **kwargs):
        """Creates a new Map instance"""
        
        maps = []
        data = []
        headers = []
    
        # convert input to maps
        for item in args:
            if isinstance(item, Map):
                maps.append(item)
            else:
                maps.append(Map.read(item))

        # sort data
        sortby = kwargs.get("sortby", "date")
        if hasattr(cls, '_sort_by_%s' % sortby):
            maps.sort(key=getattr(cls, '_sort_by_%s' % sortby)())

        # create data cube
        for map_ in maps:
            data.append(np.array(map_))
            headers.append(map_._original_header)

        obj = np.asarray(data).view(cls)
        obj._headers = headers

        return obj
    
    #pylint: disable=W0613,E1101
    def __init__(self, *args, **kwargs):
        coalign = kwargs.get("coalign", False)
        derotate = kwargs.get("derotate", False)
        
        # Coalignment
        if coalign and hasattr(self, '_coalign_%s' % coalign):
            getattr(self, '_coalign_%s' % coalign)()

        if derotate:
            self._derotate()
            
    def __array_finalize__(self, obj):
        """Finishes instantiation of the new MapCube object"""
        if obj is None:
            return

        if hasattr(obj, '_headers'):
            self._headers = obj._headers
        
    def __array_wrap__(self, out_arr, context=None):
        """Returns a wrapped instance of a MapCube object"""
        return np.ndarray.__array_wrap__(self, out_arr, context)
    
    def __getitem__(self, key):
        """Overiding indexing operation"""
        if self.ndim is 3 and isinstance(key, int):
            data = np.ndarray.__getitem__(self, key)
            header = self._headers[key]
            for cls in Map.__subclasses__():
                if cls.is_datasource_for(header):
                    return cls(data, header)

        else:
            return np.ndarray.__getitem__(self, key)
        
    def std(self, *args, **kwargs):
        """overide np.ndarray.std()"""
        return np.array(self, copy=False, subok=False).std(*args, **kwargs)
        
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
    
    def plot(self, annotate=True, axes=None, gamma=None, controls=True,
             interval=200, resample=[0.25,0.25], colorbar=False,
             **imshow_args):
        """A plot method that animates! """
        
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
        kwargs.update(imshow_args)
        
        im = axes.imshow(self[0], **kwargs)
        
        #Set current image (makes colorbar work)
        plt.sci(im)
        
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.2)
        cbar = plt.colorbar(im,cax)
        
        if resample:
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
                                            fargs = [im,annotate,ani_data],
                                            interval=interval,
                                            blit=False)
        if controls:
            axes, bax1, bax2 = plotting.add_controls(axes=axes)

            bax1._button.on_clicked(ani._start)
            bax2._button.on_clicked(ani._stop)
        
        return ani
