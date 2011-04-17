"""A Python MapCube Object

Authors: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
__author__ = "Keith Hughitt and Steven Christe"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
import pyfits

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
    def __new__(cls, input_):
        """Creates a new Map instance"""
        if isinstance(input_, str):
            headers, data = self._parseFiles(input_)

            obj = np.asarray(data).view(cls)
            obj.headers = headers
            
        elif isinstance(input_, list):
            obj = np.asarray(input_).view(cls)
        elif isinstance(input_, np.ndarray):
            obj = input_
            
        return obj
        
    def _parseFiles(self, directory):
        """Parses files in the specified directory
        
        Reads in the files at the specified location, stores their headers,
        and creates a 3d array from their contents
        """
        data    = []
        headers = []
        
        for input_ in os.listdir(directory):
            try:
                fits = pyfits.open(input_)
                data.append(fits[0].data)
                headers.append(fits[0].header)
            except IOError:
                sys.exit("Unable to read the file %s" % input_)
                
        return headers, data
            
    def __array_finalize__(self, obj):
        """Finishes instantiation of the new map object"""
        if obj is None: return

