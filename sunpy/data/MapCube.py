"""A Python MapCube Object"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import sys
import os
import pyfits
import numpy as np
from sunpy.data.sources import *
from sunpy.data.BaseMap import BaseMap
from sunpy.data.BaseMap import UnrecognizedDataSouceError

#
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
    input_ : directory, data array
        The data source used to create the map object. This can be either a
        filepath to a directory containing the images you wish to include, a 2d 
        list, or an ndarray.
    coalign : [ None | 'diff' ]
        Method to use for coalignment. If None, no coalignment will be done.
    derotate : bool
        Whether the data should be derotated

    Attributes
    ----------
    slices : list
        a list of :class:`sunpy.data.MapSlice` objects corresponding to the
        images that were used to build the MapCube.

    See Also:
    ---------
    numpy.ndarray Parent class for the MapCube object
    :class:`sunpy.Map`
        
    Examples
    --------
    >>> mapcube = sunpy.MapCube('images/')
    >>> mapcube.slices[0].header['crpix1']
    2050.6599120000001
    """
    def __new__(cls, input_, *args):
        """Creates a new Map instance"""
        if isinstance(input_, str):
            data = []
            slices = []
            
            for filename in os.listdir(input_):
                try:
                    fits = pyfits.open(os.path.join(input_, filename))
                    data.append(fits[0].data)
                    slices.append(cls.parse_header(fits[0].header))
                except IOError:
                    sys.exit("Unable to read the file %s" % filename)

            obj = np.asarray(data).view(cls)
            obj.slices = slices
            
        elif isinstance(input_, list):
            obj = np.asarray(input_).view(cls)
        elif isinstance(input_, np.ndarray):
            obj = input_

        return obj
    
    def __init__(self, input_, coalign=False, derotate=False):
        # Coalignment
        if coalign and hasattr(self, '_coalign_%s' % coalign):
            getattr(self, '_coalign_%s' % coalign)()

        if derotate:
            obj._derotate()
        
    
    @classmethod
    def parse_header(cls, header):
        """Returns a MapSlice instance corresponding to an image header.
        
        Attempts to construct a `MapSlice` instance using the header information
        provided. `MapSlice` instances (e.g. `AIAMapSlice`) are basically just
        empty (data-less) Map objects. This provides a way to keep track of
        the meta-information for each of the images that were used to build
        the MapCube separately from the data.
        
        Parameters
        ----------
        header : dict
            The image header for which a `MapSlice` should be built
        """
        for cls in BaseMap.__subclasses__():
            if cls.is_datasource_for(header):
                return cls.as_slice(header)
        raise UnrecognizedDataSouceError
    
    def _derotate(self):
        """Derotates the layers in the MapCube"""
        pass
            
    def __array_finalize__(self, obj):
        """Finishes instantiation of the new map object"""
        if obj is None: return
        
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
    
