"""A Python MapCube Object

Authors: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
import pyfits
from BaseMap import BaseMap
from BaseMap import UnrecognizedDataSouce
from sources import *

#
# 2011/04/13: Should Map be broken up into Map and MapHeader classes? This way
# mapped header values can be used in MapCube without having to keep extra
# copies of the data..
#
class MapCube(np.ndarray):
    """
    MapCube(input)
    
    A spatially-aware data array based on the SolarSoft Map object

    Attributes
    ----------
    headers : list
        a list of dictionaries for each image header in the collection

    Parameters
    ----------
    input_ : directory, data array
        The data source used to create the map object. This can be either a
        filepath to a directory containing the images you wish to include, a 2d 
        list, or an ndarray.
        
    See Also:
    ---------
    numpy.ndarray Parent class for the MapCube object
    :class:`sunpy.Map`
        
    Examples
    --------
    >>> mapcube = sunpy.MapCube('images/')

    """
    def __new__(cls, input_, coalign=False, derotate=False):
        """Creates a new Map instance"""
        if isinstance(input_, str):
            headers, data = self._parseFiles(input_)

            obj = np.asarray(data).view(cls)
            obj.slices = slices
            
        elif isinstance(input_, list):
            obj = np.asarray(input_).view(cls)
        elif isinstance(input_, np.ndarray):
            obj = input_
            
        if coalign:
            self._coalign()
        if derotate:
            self.derotate()
            
        return obj
        
    def _parseFiles(self, directory):
        """Parses files in the specified directory
        
        Reads in the files at the specified location, stores their headers,
        and creates a 3d array from their contents
        """
        data = []
        slices = []
        
        for input_ in os.listdir(directory):
            try:
                fits = pyfits.open(input_)
                data.append(fits[0].data)
                slices.append(self._parse_header(fits[0].header))                
            except IOError:
                sys.exit("Unable to read the file %s" % input_)
                
        return slices, data
    
    def _parse_header(self):
        """Returns a MapSlice instance corresponding to an image header.
        
        Attempts to construct a MapSlice instance using the header information
        provided. MapSlice instances (e.g. AIAMapSlice) are basically just
        empty (data-less) Map objects. This provides a way to keep track of
        the meta-information for each of the images that were used to build
        the MapCube separately from the data. 
        """
        for cls in BaseMap.__subclasses__():
            if cls.is_datasource_for(fits[0].header):
                return cls.as_slice(header)
        raise UnrecognizedDataSouce
    
    def _coalign(self):
        """Coaligns the layers in the MapCube"""
        pass
    
    def derotate(self):
        """Derotates the layers in the MapCube"""
        pass
            
    def __array_finalize__(self, obj):
        """Finishes instantiation of the new map object"""
        if obj is None: return
