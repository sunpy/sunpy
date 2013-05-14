"""A Python MapCube Object"""
from __future__ import absolute_import
#pylint: disable=W0401,W0614,W0201,W0212,W0404

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import copy
from datetime import timedelta

from sunpy.map import Map
from sunpy.map.sources import *
from sunpy.lightcurve import LightCurve
from sunpy.util import plotting
from sunpy.coords import rot_hpc

__all__ = ['MapCube']

# (https://github.com/sunpy/sunpy/issues/397)
# 2011/04/13: Should Map be broken up into Map and MapHeader classes? This way
# mapped header values can be used in MapCube without having to keep extra
# copies of the data..
# 2013/04/14
# The mapcube has to be easily accessible along all dimensions.  At the moment.
# it is not.
# MC = MapCube(list of fits files, array)
# 
# 1-d (lightcurve)
# lc = MC.sub(0, 0, 100:200)
#
# 2-d (map)
# mp = MC.sub(100:200, 200:300, 4)
# 
# 3-4 (mapcube)
# mc = MC.sub(100:200, 200:300, 4:40)
#
# Perhaps NDData is the answer?
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

    ordering : {"order":  # one-dimensional array of orderable objects,
                "description": # a description of the meaning of "order"
                "units": } #the units of "order"
    Attributes
    ----------
    headers : list
        a list of dictionaries containing the original and normalized header
        tags for the files used to build the MapCube.

    ordering : dictionary
        a dictionary that contains the following tags
            "order": the numerical list that orders the maps in the mapcube
            "description": a description of the meaning of "order"
            "units": the units of "order"

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

        # default ordering
        ordering = kwargs.get("ordering", {"order": range(0, len(maps)),
                                          "description": "default ordering",
                                          "units": 'index'})

        # sort data.  a sort method overwrites the existing ordering
        sortby = kwargs.get("sortby", "date")
        if hasattr(cls, '_sort_by_%s' % sortby):
            sort_key = getattr(cls, '_sort_by_%s' % sortby)()
            maps.sort(key=sort_key)
            ordering = {"order": [sort_key(map_) for map_ in maps],
                        "description": 'time',
                        "units": ''}

        # create data cube
        for map_ in maps:
            data.append(np.array(map_))
            headers.append(map_._original_header)

        obj = np.asarray(data).view(cls)
        obj._headers = headers
        obj._ordering = ordering

        return obj

    #pylint: disable=W0613,E1101
    def __init__(self, *args, **kwargs):
        coalign = kwargs.get("coalign", False)
        derotate = kwargs.get("derotate", False)

        # coalignment
        if coalign and hasattr(self, '_coalign_%s' % coalign):
            getattr(self, '_coalign_%s' % coalign)()

        if derotate is not False:
            if hasattr(self, '_derorate_%s' % derotate):
                getattr(self, '_derotate_%s' % derotate)(**kwargs)

    def __array_finalize__(self, obj):
        """Finishes instantiation of the new MapCube object"""
        if obj is None:
            return

        if hasattr(obj, '_headers'):
            self._headers = obj._headers

        if hasattr(obj, '_ordering'):
            self._ordering = obj._ordering

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

    def get_lightcurve_by_array_index(self, x, y):
        """Returns a lightcurve object at a given pixel"""
        data = [m[x, y] for m in self]
        return LightCurve.create({m.name: data}, index=self.ordering["order"])

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
        return lambda m: m.date  # maps.sort(key=attrgetter('date'))

    def _derotate_by_latitude(self, index=0, use_order=False):
        """Derotates the layers in the MapCube.  Derotates each image using
        the latitudinal dependence defined by diff_rot.  Derotates the stack of
        images to the map in the stack at position 'index'.  If use_order is
        True then we assume that index is of the same type as ordering["order"]
        and the map stack is derotated to the closest map in the stack to the
        value of index passed.  Another way of putting this is to say 'for the
        mapcube find the map in mapcube that minimizes
        abs(map._ordering["order"] - index).  If use_order is False then the
        maps in the mapcube are derotated relative to mapcube[i]."""
        pass


    def derotate_by_center_of_fov(self, **kwargs):
        """Derotate layers of the MapCube using the center of the FOV in each
        layer only. Should be faster than _derotate_by_latitude.   Derotates
        the stack of images to the map in the stack at position 'index'.
        If use_order = True then we assume that index is of the same type as
        ordering["order"] and the map stack is derotated to the closest map in
        the stack to the value of index passed.  Another way of putting this is
        to say 'for the mapcube, find the map in mapcube that minimizes
        abs(map._ordering["order"] - index).  If use_order is False then the
        maps in the mapcube are derotated relative to mapcube[i]."""

        index = kwargs.get("index", 0)
        use_order = kwargs.get("use_order", False)
        if use_order:
            difference = index - self._ordering["order"]
            if isinstance(difference, timedelta):
                difference = np.absolute((index - self._ordering["order"]).
                                         to_seconds())
            else:
                difference = np.absolute((index - self._ordering["order"]))
            index = np.where(difference == difference.min())[0][0]
        print index

        # Center of the field of view of the base map
        xcen_base = self._headers[index]["xcen"]
        ycen_base = self._headers[index]["ycen"]

        # Pixel size of the pixels in the base map
        xpixel_size = self._headers[index]["CDELT1"]
        ypixel_size = self._headers[index]["CDELT2"]

        # for each map calculate the derotation value in pixels, move the data,
        # recreate a map with the moved data, and create a MapCube.
        for m in self:
            new_xcen, new_ycen = rot_hpc(m, tstart=m.header["date_obs"],
                                 interval = self[index].header["date_obs"])
            xpixel_difference = (new_xcen - xcen_base)/xpixel_size
            ypixel_difference = (new_ycen - ycen_base)/ypixel_size
            
            pass

        return

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
