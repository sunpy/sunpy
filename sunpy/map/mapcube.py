"""A Python MapCube Object"""
#pylint: disable=W0401,W0614,W0201,W0212,W0404

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap, UnrecognizedDataSouceError
from sunpy.map.sources import * #@UnusedWildImport
import numpy as np
import os
import pyfits

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
    input_ : directory, glob string, data array, 
        The data source used to create the map object. This can be either a
        filepath to a directory containing the images you wish to include, a
        globbable string such as "data/*.fits", a 2d list, or an ndarray.
    coalign : [ None | 'diff' ]
        Method to use for coalignment. If None, no coalignment will be done.
    derotate : bool
        Whether the data should be derotated

    Attributes
    ----------
    headers : list
        a list of dictionaries containing the original and normalized header
        tags for the files used to build the MapCube.

    See Also
    --------
    numpy.ndarray Parent class for the MapCube object
    :class:`sunpy.map.BaseMap`
        
    Examples
    --------
    >>> mapcube = sunpy.MapCube('images/')
    >>> mapcube[0].plot()
    >>> mapcube[3].header.get('crpix1')
    2050.6599120000001
    """
    def __new__(cls, input_, sortby="date"):
        """Creates a new Map instance"""
        
        # Directory of files
        if isinstance(input_, basestring):
            filepaths = []
            fits_arr = []
            data = []
            headers = []

            # directory
            if os.path.isdir(input_):
                for filename in os.listdir(input_):
                    filepaths.append(os.path.join(input_, filename))

            # glob string
            else:
                from glob import glob
                filepaths = glob(input_)
                
            # read in files
            for filepath in filepaths:
                fits = pyfits.open(filepath)
                
                # append normalized header tags for use during sorting
                found_header_match = False
                
                for subcls in BaseMap.__subclasses__(): #pylint: disable=E1101
                    if subcls.is_datasource_for(fits[0].header):
                        found_header_match = True
                        fits.norm_header = subcls.get_properties(fits[0].header)
                if not found_header_match:
                    raise UnrecognizedDataSouceError

                fits_arr.append(fits)

            # sort data
            if sortby and hasattr(cls, '_sort_by_%s' % sortby):
                fits_arr.sort(key=getattr(cls, '_sort_by_%s' % sortby)())

            # create data cube
            for fits in fits_arr:
                data.append(fits[0].data)
                headers.append(fits[0].header)

            obj = np.asarray(data).view(cls)
            obj._headers = headers

        # List of data or filepaths
        elif isinstance(input_, list):
            obj = np.asarray(input_).view(cls)

        # ndarray
        elif isinstance(input_, np.ndarray):
            obj = input_

        return obj
    
    #pylint: disable=W0613,E1101
    def __init__(self, input_, coalign=False, derotate=False):
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
            for cls in BaseMap.__subclasses__():
                if cls.is_datasource_for(header):
                    return cls(data, header)
            raise UnrecognizedDataSouceError
        else:
            return np.ndarray.__getitem__(self, key)
        
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
        #from operator import attrgetter
        #return attrgetter("norm_header.date")
        return lambda fits: fits.norm_header.get("date")
    
    def _derotate(self):
        """Derotates the layers in the MapCube"""
        pass
    
    def plot(self):
        """A basic plot method (not yet implemented)"""
        pass
    
if __name__ == "__main__":
    import sunpy
    m = sunpy.MapCube("/home/hughitt1/Dropbox/eitwave")
    #print(m[0,0:5,0:5])
    #print()
    repr(m[2].base)
