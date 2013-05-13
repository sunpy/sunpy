from __future__ import absolute_import

import os
import glob
import urllib2

import numpy as np

from sunpy.map.map import Map
from sunpy.map.header import MapHeader
from sunpy.map.mapcube import MapCube
from sunpy.map.compositemap import CompositeMap
from sunpy.map.sources import *

from sunpy.util.datatype_factory_base import RegisteredFactoryBase

__all__ = ['make_map', 'read_header']


MapFactoryArgParser(object):

    # This can be hacked later to use ConditionalDispatch, but I'm not going to
    # do that.  (This explains why it is a class and not a function at the moment.  Really, this should probably be a class method of Map (the factory).

    @classmethod
    def read_file(cls, fname):
        
        # File gets read here.  This needs to be generic enough to seamlessly call a fits file or a jpeg2k file, etc
        
        data = None
        meta = None
        
        return (data, header)

    @classmethod
    def __call__(cls, *args, **kwargs):
        
        data_header_pairs = list()
        already_maps = list()
        
        # For each of the arguments, handle each of the cases
        i = 0
        while i < len(args):
            
            arg = args[i]
            
            # Data-header pair in a tuple
            if (type(arg) in [tuple, list]) and 
                          isinstance(arg[0],np.ndarray) and # or NDData or something else?
                          isinstance(arg[1],SunpyMetaBase): # FITSHeader, JP2kHeader, OrderedDict, dict?
                data_header_pairs.append(arg)
            
            # Data-header pair not in a tuple
            elif isinstance(arg, np.ndarray) and # or NDData or something else?
                          isinstance(args[i+1],SunpyMetaBase): # FITSHeader, JP2kHeader, OrderedDict, dict? 
                pair = (args[i], args[i+1])
                data_header_pairs.append(pair)
                i += 1 # an extra increment to account for the data-header pairing
            
            # File name
            elif type(arg) is basestring and 
                          os.path.isfile(os.path.expanduser(arg)):
                path = os.path.expanduser(arg)
                pair = cls.read_file(path)
                data_header_pairs.append(pair)
            
            # Directory
            elif type(arg) is basestring and 
                          os.path.isdir(os.path.expanduser(arg)):
                path = os.path.expanduser(arg)
                files = [os.path.join(directory, elem) for elem in os.listdir(path)]
                data_header_pairs += map(cls.read_file, files)
            
            # Glob
            elif type(arg) is basestring and 
                          '*' in arg:
                files = glob.glob( os.path.expanduser(arg) )
                data_header_pairs += map(cls.read_file, files)
            
            # Already a Map
            elif isinstance(arg, MapBase):
                already_maps.append(arg)
                
            # A URL
            elif type(arg) is basestring and 
                          urllib2.urlopen(arg):
                default_dir = sunpy.config.get("downloads", "download_dir")
                path = download_file(url, default_dir)
                pair = cls.read_file(path)
                data_header_pairs.append(pair)
        
            i += 1
            
        return data_header_pairs, already_maps

class Map(RegisteredFactoryBase):
    
    GenericWidgetType = GenericMap
    
    def __new__(cls, composite=False, cube=False, *args, **kwargs):

        if cls is Map:
            
            # Get list of data-header pairs, e.g., [(d1, h1), (d2, h2), ...]
            data_header_pairs = MapFactoryArgParser(*args, **kwargs)
            
            # If the list is meant to be a cube, instantiate a map cube
            if cube:
                return MapCube(*data_header_pairs, **kwargs)

            # If the list is meant to be a composite mape, instantiate one
            if composite:
                return CompositeMap(*data_header_pairs, **kwargs)
                
            # Otherwise, each pair in the list gets built on its own
            new_maps = list()
            
            for data, header in zip(*data_header_pairs):
                
                # Test to see which type of Map this pair is.  If none of the
                # registered Map types match, use a generic map.
                WidgetType = None
                for key in self.registry:
                    
                    if cls.registry[key](data, header, **kwargs):
                        WidgetType = key
                        break
                else:
                    WidgetType = cls.GenericWidgetType
                    
                # Instantiate the new map.
                new_maps.append(WidgetType(data, header, **kwargs))
                
            # If there was only one map instantiated, return that, otherwise
            # return the list of them.
            return new_maps[1] if len(new_maps) == 1 else new_maps
        else:
            return super(Map, cls).__new__(cls, *args, **kwargs)
 


def make_map(*args, **kwargs):
    """Processes one or more inputs and returns a Map, MapCube, or CompositeMap
    instance.

    Parameters
    ----------
    args : filepath(s), data array
        The data source used to create the map object. This can be either a
        filepath to an image, a 2d list, or an ndarray.
    type : {'composite' | 'cube'}
        Type of multimap to construct when passed more than one input. The
        default choice is a CompositeMap which is more lenient with respect
        to how similar the input data is.

    Returns
    -------
    out : Map, MapCube, CompositeMap
        Returns a  subclass instance

    Examples
    --------
    >>> import sunpy
    >>> sunpy.make_map("file.fts")
    >>> sunpy.make_map("file1.fts", "file2.fts",..)
    >>> sunpy.make_map(["file1.fts", "file2.fts",..])
    >>> sunpy.make_map("path/to/files/*.fts")
    >>> sunpy.make_map(Map)
    >>> sunpy.make_map(Map1, Map2,..)
    >>> sunpy.make_map([[0, 1],[2, 3]], {'telescop': 'sunpy',..})

    """
    if len(args) is 0:
        raise TypeError("Invalid input.")

    # First check to see if data/header were passed in    
    if isinstance(args[0], list) or isinstance(args[0], np.ndarray):
        data = None

        # n-dimensional list
        if isinstance(args[0][0], list) or isinstance(args[0], np.ndarray):
            data = args[0]
        else:
            try:
                float(args[0][0])
            except (ValueError, TypeError):
                pass
            else:
                # 1-dimensional data
                data = args[0]

        # if either of the above cases hold, then create a new Map
        if data is not None:
            if len(args) > 1:
                return Map(args[0], args[1])
            else:
                return Map(args[0], {})


    # If not, check for one or more maps or filepaths
    if len(args) == 1:
        # String
        if isinstance(args[0], basestring):
            filepath = os.path.expanduser(args[0])

            # Wildcard string
            if filepath.find("*") != -1:
                maps = glob.glob(filepath)
            # Directory (use all files)
            elif os.path.isdir(filepath):
                maps = [os.path.join(filepath, x) for x in os.listdir(filepath)]

            # Filepath
            else:
                return Map.read(filepath)

        # Map/MapCube/CompositeMap
        elif (isinstance(args[0], Map) or 
              isinstance(args[0], CompositeMap) or 
              isinstance(args[0], MapCube)):
            return args[0]

        # List of filepaths or Maps
        elif isinstance(args[0], list):
            # list of maps or filepaths
            maps = args[0]

        # Unrecognized input
        else:
            raise InvalidMapInput("Invalid input for make_map. Please specify "
                                  "one or more filepaths, Maps, directories, "
                                  "or wildcard expressions.")
    else:
        maps = args

    # Make sure we found some data
    if len(maps) is 0:
        raise NoMapsFound("Specified path contains no valid files.")

    mtype = kwargs.get("type", "composite")

    # MapCube
    if mtype == "cube":
        return MapCube(*maps)
    # CompositeMap (default)
    elif mtype == "composite":
        return CompositeMap(*maps)
    else:
        raise InvalidMapType("Invalid multi-map type specified. Please choose "
                             "between 'composite' or 'cube'.")

def read_header(filepath):
    """Parses a file header and return some important parameters"""
    return Map.read_header(filepath)

class InvalidMapInput(ValueError):
    """Exception to raise when input variable is not a Map instance and does
    not point to a valid Map input file."""
    pass

class InvalidMapType(ValueError):
    """Exception to raise when an invalid type of map is requested with make_map
    """
    pass

class NoMapsFound(ValueError):
    """Exception to raise when input does not point to any valid maps or files 
    """
    pass